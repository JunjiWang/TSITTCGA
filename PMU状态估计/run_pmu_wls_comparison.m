function summary=run_pmu_wls_comparison()
%RUN_PMU_WLS_COMPARISON Paired noise-only / sparse-bad-data accuracy study.
% All 45 frozen cases of each IEEE14/118/300 system are retained. Parameters,
% truth and Gaussian noise are fixed, with no tuning, restarts or new noise.
% WLS is OLS because the supplied covariance is sigma^2 I. Its LS validity
% is checked separately; the robust CNE criterion is NOT applied to WLS.
% Each robust method is freshly solved once per condition. This deterministic
% accuracy study does not replace the existing repeated runtime comparison.
here=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;
threads=maxNumCompThreads(1);cleanup=onCleanup(@()restore(oldpath,oldrng,threads)); %#ok<NASGU>
addpath(here,'-begin');c=pmu_parameters();
helpers={'pmu_parameters','pmu_run_cgpm','pmu_prepare_backend', ...
    'pmu_map','huber_irls_same_model','TSITTCGA_precision'};
for h=helpers
    assert(strcmpi(fileparts(which(h{1})),here),'PMU:Shadowing','Wrong source: %s',h{1});
end
out=fullfile(here,'results');if ~isfolder(out),mkdir(out);end
prefix=['wls_paired_' c.profile];matfile=fullfile(out,[prefix '.mat']);
assert(~isfile(matfile),'PMU:ExistingResults','Preserve existing results before rerunning: %s',matfile);
protocol=struct('Created',char(datetime('now')),'MATLAB',version,'Config',c, ...
    'Systems',{c.systems},'Threads',1,'CasesPerSystem',45,'Completed',false, ...
    'Scenarios',{{'GaussianOnly','SparseBadData'}}, ...
    'Pairing','Same G, true state and stored noise; add the original stored bad vector.', ...
    'WLS','sparse(G)\y; covariance sigma^2 I, no robust-equation stopping test.', ...
    'AccuracyOnly',true,'TimingComparison',false,'NoTuning',true,'FreshRunsPerCondition',1, ...
    'RobustTolerance',c.target,'RobustBudget',c.budget,'InputFiles',struct([]));
files=dir(fullfile(here,'*.m'));snapshot=struct([]);
for j=1:numel(files)
    snapshot(j).Name=files(j).name;snapshot(j).Text=fileread(fullfile(here,files(j).name)); %#ok<AGROW>
end
protocol.SourceSnapshot=snapshot;
rows=struct([]);endpoints={};paired=struct([]);completed=false;
save(matfile,'protocol','rows','endpoints','paired','completed','-v7.3');
for sy=1:numel(c.systems)
    sys=c.systems{sy};inputfile=fullfile(here,'data',['input_' sys '.mat']);
    d=load(inputfile,'instances');cases=d.instances;
    assert(numel(cases)==45,'PMU:Coverage','Expected all 45 cases.');
    protocol.InputFiles(sy).System=sys;protocol.InputFiles(sy).Path=inputfile;
    protocol.InputFiles(sy).SHA256=file_sha256(inputfile);
    for k=1:numel(cases)
        raw=cases{k};clean_y=raw.ytrue(:)+raw.noise(:);
        assert(norm(raw.y(:)-clean_y-raw.bad(:),inf)<=1e-12*max(1,norm(raw.y,inf)), ...
            'PMU:Pairing','Original y does not equal truth + stored noise + bad data.');
        assert(norm(raw.G*raw.true_xi(:)-raw.ytrue(:))<=1e-10*max(1,norm(raw.ytrue)), ...
            'PMU:Truth','Measurement truth does not match G*true_xi.');
        assert(raw.sigma>0&&all(isfinite(raw.noise)),'PMU:Noise','Invalid stored Gaussian noise.');
        rr=cell(2,3);ee=nan(2,3);
        for sc=1:2
            inst=raw;if sc==1,inst.y=clean_y;end
            backend=c.backend;
            if isfield(c,'backend_by_system'),backend=c.backend_by_system.(sys);end
            [inst,~]=pmu_prepare_backend(inst,backend);
            G=sparse(inst.G);xi_wls=G\inst.y;fit=inst.y-G*xi_wls;
            % Data-scale normalization includes rounding in y-G*xi. Scaling
            % only by ||fit|| spuriously rejects an accurate high-scale LS.
            normal=norm(G'*fit,inf)/max(1,norm(G,1)*(norm(inst.y,inf)+norm(G,inf)*norm(xi_wls,inf)));
            projection_error=norm(G*xi_wls-raw.Q*(raw.Q'*inst.y))/max(1,norm(inst.y));
            assert(all(isfinite(xi_wls))&&normal<=1e-12&&projection_error<=1e-11, ...
                'PMU:WLS','WLS data-scaled or independent projection validation failed.');
            rr{sc,1}=struct('xi',xi_wls,'relerr',norm(xi_wls-inst.true_xi)/norm(inst.true_xi), ...
                'success',true,'status','least_squares_solution','residual',NaN, ...
                'feasibility',NaN,'obj',NaN,'Itr',NaN,'NF',NaN,'LSNormalResidualScaled',normal, ...
                'LSProjectionRelativeError',projection_error);
            for a=2:3
                if a==2
                    r=pmu_run_cgpm(inst,c.solvers{1},c.params{1},c.lambda,c.initial);
                else
                    r=huber_irls_same_model(inst,c.lambda,c.target,c.budget,false);
                end
                r=independent_check(inst,r,c,a==2);rr{sc,a}=r;
            end
            if sc==1,scenario='GaussianOnly';else,scenario='SparseBadData';end
            names={'WLS','TSITTCGA','Huber-IRLS'};
            for a=1:3
                r=rr{sc,a};ee(sc,a)=r.relerr;
                if a==1,reached=NaN;ls=r.LSNormalResidualScaled;pe=r.LSProjectionRelativeError;
                else,reached=double(r.success);ls=NaN;pe=NaN;end
                row=struct('System',sys,'Case',k,'Scenario',scenario,'Method',names{a}, ...
                    'Measurements',inst.m,'StateDimension',size(G,2),'CNEDimension',2*inst.m, ...
                    'Currents',raw.Currents,'BadDataIfContaminated',nnz(raw.bad),'SeedIndex',raw.SeedIndex, ...
                    'RandomSeed',raw.RandomSeed,'ValidState',all(isfinite(r.xi)), ...
                    'RobustSolved',reached,'Status',r.status,'RelativeStateError',r.relerr, ...
                    'RelativeStateErrorPercent',100*r.relerr,'OriginalResidual',r.residual, ...
                    'Feasibility',r.feasibility,'Objective',r.obj,'LSNormalResidualScaled',ls, ...
                    'LSProjectionRelativeError',pe);
                if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
                endpoints{end+1}=r; %#ok<AGROW>
            end
        end
        allrobust=all(cellfun(@(r)r.success,rr(:,2:3)),'all');
        pr=struct('System',sys,'Case',k,'Currents',raw.Currents,'BadData',nnz(raw.bad), ...
            'SeedIndex',raw.SeedIndex,'RandomSeed',raw.RandomSeed,'CommonRobustSuccess',allrobust, ...
            'WLSCleanRE',ee(1,1),'WLSBadRE',ee(2,1),'TSCleanRE',ee(1,2),'TSBadRE',ee(2,2), ...
            'IRLSCleanRE',ee(1,3),'IRLSBadRE',ee(2,3), ...
            'TSBetterThanWLSBad',ee(2,2)<ee(2,1),'IRLSBetterThanWLSBad',ee(2,3)<ee(2,1), ...
            'TSBadToWLSBadRatio',ee(2,2)/max(ee(2,1),realmin));
        if isempty(paired),paired=pr;else,paired(end+1)=pr;end %#ok<AGROW>
        fprintf('%s case %2d/45 WLS clean/bad %.4g/%.4g TS %.4g/%.4g IRLS %.4g/%.4g common=%d\n', ...
            sys,k,ee(1,1),ee(2,1),ee(1,2),ee(2,2),ee(1,3),ee(2,3),allrobust);
        save(matfile,'protocol','rows','endpoints','paired','completed','-v7.3');
        writetable(struct2table(rows),fullfile(out,[prefix '_runs.csv']));
        writetable(struct2table(paired),fullfile(out,[prefix '_paired.csv']));
    end
end
completed=true;protocol.Completed=true;
summary=make_summary(rows,paired,c.systems);
strata=table();
for sy=1:numel(c.systems)
    pr=paired(strcmp({paired.System},c.systems{sy}));
    groups=unique([[pr.Currents].',[pr.BadData].'],'rows');
    for j=1:size(groups,1)
        ids=[pr([pr.Currents]==groups(j,1)&[pr.BadData]==groups(j,2)).Case];
        subrows=rows(strcmp({rows.System},c.systems{sy})&ismember([rows.Case],ids));
        subpair=pr(ismember([pr.Case],ids));
        st=make_summary(subrows,subpair,c.systems(sy));
        st.Currents=repmat(groups(j,1),height(st),1);st.BadData=repmat(groups(j,2),height(st),1);
        strata=[strata;st]; %#ok<AGROW>
    end
end
save(matfile,'protocol','rows','endpoints','paired','summary','strata','completed','-v7.3');
writetable(strata,fullfile(out,[prefix '_strata.csv']));
writetable(summary,fullfile(out,[prefix '_summary.csv']));disp(summary);
end

function r=independent_check(inst,r,c,is_ts)
m=inst.m;s=r.s(:);assert(numel(s)==2*m,'Invalid split endpoint.');
if ~all(isfinite(s))
    assert(~r.success,'PMU:NonfiniteSuccess','Nonfinite endpoint claimed success.');
    r.IndependentResidual=Inf;r.IndependentObjective=NaN;r.IndependentRelativeError=NaN;
    r.IndependentStateDifference=NaN;r.IndependentVerified=false;return;
end
e=s(1:m)-s(m+1:end);assert(norm(e-r.e)<=1e-12*max(1,norm(e)),'Split mismatch.');
xi=sparse(inst.G)\(inst.y-e);fit=inst.y-e-inst.G*xi;
fv=[min(s(1:m),c.lambda-fit);min(s(m+1:end),c.lambda+fit)];
fn=norm(fv);feas=norm(min(s,0));obj=.5*dot(fit,fit)+c.lambda*norm(e,1);
re=norm(xi-inst.true_xi)/norm(inst.true_xi);
if ~all(isfinite([fn;feas;obj;re;xi]))
    assert(~r.success,'PMU:NonfiniteSuccess','Nonfinite recomputation claimed success.');
    r.IndependentResidual=fn;r.IndependentObjective=obj;r.IndependentRelativeError=re;
    r.IndependentStateDifference=NaN;r.IndependentVerified=false;return;
end
assert(abs(fn-r.residual)<=1e-9&&abs(obj-r.obj)<=1e-10*max(1,abs(obj)) ...
    &&abs(re-r.relerr)<=1e-10*max(1,re),'PMU:IndependentEndpoint','Fresh LS check disagrees.');
if r.success
    assert(fn<=c.target&&feas<=c.feasibility_target&&all(isfinite(xi)), ...
        'PMU:IndependentSuccess','Successful endpoint fails the original tolerance.');
end
if is_ts,assert(r.count_ok&&r.report_ok,'PMU:NF','NF or residual reporting mismatch.');
else,assert(r.monotonicity_violations==0,'PMU:IRLS','IRLS objective increased.');end
r.IndependentResidual=fn;r.IndependentObjective=obj;r.IndependentRelativeError=re;
r.IndependentStateDifference=norm(xi-r.xi);r.IndependentVerified=true;
end

function T=make_summary(rows,paired,systems)
out=struct([]);names={'WLS','TSITTCGA','Huber-IRLS'};scenarios={'GaussianOnly','SparseBadData'};
for sy=1:numel(systems)
    common=[paired(strcmp({paired.System},systems{sy})).CommonRobustSuccess];
    ids=[paired(strcmp({paired.System},systems{sy})).Case];commonids=ids(common);
    for sc=1:2
        for a=1:3
            z=rows(strcmp({rows.System},systems{sy})&strcmp({rows.Scenario},scenarios{sc})&strcmp({rows.Method},names{a}));
            v=[z.RelativeStateErrorPercent];v=v(isfinite(v));
            cv=[z(ismember([z.Case],commonids)).RelativeStateErrorPercent];
            if a==1,solved=NaN;else,solved=sum([z.RobustSolved]);end
            q=prctile(v,[25 50 75]);
            r=struct('System',systems{sy},'Scenario',scenarios{sc},'Method',names{a}, ...
                'Cases',numel(z),'ValidStates',sum([z.ValidState]),'RobustSolved',solved, ...
                'MeanREPercentAll',mean(v),'MedianREPercentAll',q(2),'Q1REPercentAll',q(1),'Q3REPercentAll',q(3), ...
                'CommonBothScenariosCases',numel(commonids),'MeanREPercentCommon',mean(cv));
            if isempty(out),out=r;else,out(end+1)=r;end %#ok<AGROW>
        end
    end
end
T=struct2table(out);
end

function s=file_sha256(p)
digest=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(p,'rb');
assert(fid>=0,'Cannot read input for hashing.');cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid),b=fread(fid,1048576,'*uint8');digest.update(typecast(b,'int8'));end
s=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function restore(p,r,t)
path(p);rng(r);maxNumCompThreads(t);
end
