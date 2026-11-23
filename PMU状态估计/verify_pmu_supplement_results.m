function report=verify_pmu_supplement_results()
%VERIFY_PMU_SUPPLEMENT_RESULTS Fresh LS audit of the two supplemental studies.
% No solver, saved projector, prepared recovery object or timing ranking is
% used here. All states are recovered with fresh sparse G\ right-hand sides.
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');
oldpath=path;threads=maxNumCompThreads(1);
cleanup=onCleanup(@()restore_session(oldpath,threads)); %#ok<NASGU>
addpath(here,'-begin');c=pmu_parameters();
wfile=fullfile(out,['wls_paired_' c.profile '.mat']);
qfile=fullfile(out,['qp_validation_' c.profile '.mat']);
w=load(wfile);q=load(qfile);
assert(w.completed&&q.completed,'PMUAudit:Incomplete','Complete both studies before auditing.');
assert(numel(w.rows)==810&&numel(w.endpoints)==810&&numel(w.paired)==135, ...
    'PMUAudit:Coverage','Expected 810 WLS-study endpoints and 135 pairs.');
assert(numel(q.rows)==10&&isequal(size(q.endpoints),[10 2])&&height(q.selection)==10, ...
    'PMUAudit:QPCoverage','Expected ten QP comparisons and twenty endpoints.');
check_config(w.protocol.Config,c);check_config(q.protocol.Config,c);
assert(w.protocol.AccuracyOnly&&~w.protocol.TimingComparison ...
    &&q.protocol.QPReferenceOnly&&~q.protocol.SamePrecisionTimingComparison);
% Parameters are checked numerically above. The projector is checked by
% fresh least-squares residuals below, so later equivalent implementations
% do not invalidate archived accuracy data. Compare the numerical kernels.
core={'TSITTCGA_precision.m','pmu_run_cgpm.m', ...
    'pmu_map.m','huber_irls_same_model.m'};
source_checks=check_sources(w.protocol.SourceSnapshot,core,here);
source_checks=source_checks+check_sources(q.protocol.SourceSnapshot, ...
    [core,{'pmu_qp_reference.m'}],here);
fpath=fullfile(out,sprintf('final_budget%d_rep%d_ieee14_ieee118_%s.mat', ...
    c.budget,c.repetitions,c.profile));
formal_available=isfile(fpath);
if formal_available
    f=load(fpath,'rows','endpoints','completed','protocol');
    assert(f.completed,'PMUAudit:FormalIncomplete','The matching formal comparison is incomplete.');
    check_config(f.protocol.Config,c);
else
    fprintf('No matching historical efficiency file; auditing accuracy endpoints independently.\n');
end
stats=struct('WLSChecks',0,'RobustChecks',0,'RobustSuccesses',0, ...
    'RetainedRobustFailures',0,'NonfiniteFailedEndpoints',0,'NativeNFChecks',0, ...
    'FormalReproductionChecks',0,'FormalExactReproductionChecks',0, ...
    'FormalCrossBackendChecks',0,'PairChecks',0,'QPChecks',0,'QPVerified',0, ...
    'RetainedQPFailures',0,'QPKKTChecks',0,'QPDualityChecks',0,'InputHashChecks',0, ...
    'MaxWLSNormalDataScaled',0,'MaxFreshResidualDifference',0, ...
    'MaxFreshObjectiveDifference',0,'MaxFreshRelativeErrorDifference',0, ...
    'MaxFreshStateDifference',0,'MaxFreshSuccessfulResidual',0, ...
    'MaxFormalEndpointRelativeDifference',0,'MaxFreshQPKKTScaled',0, ...
    'MaxFreshQPScaledDualityGap',0);
checked=false(1,810);fresh_re=nan(1,810);fresh_success=false(1,810);
qp_checked=false(1,10);systems=c.systems;names={'WLS','TSITTCGA','Huber-IRLS'};
for sy=1:numel(systems)
    sys=systems{sy};input_path=fullfile(here,'data',['input_' sys '.mat']);
    input_hash=sha256_file(input_path);
    iw=find(strcmp({w.protocol.InputFiles.System},sys));
    iq=find(strcmp({q.protocol.InputFiles.Name},['input_' sys '.mat']));
    assert(isscalar(iw)&&isscalar(iq)&&strcmp(input_hash,w.protocol.InputFiles(iw).SHA256) ...
        &&strcmp(input_hash,q.protocol.InputFiles(iq).SHA256),'PMUAudit:InputHash','Input hash changed.');
    stats.InputHashChecks=stats.InputHashChecks+2;
    d=load(input_path,'instances');assert(numel(d.instances)==45);
    chosen=double(q.protocol.SelectedCaseIDs{sy}(:));
    actual=double(q.selection.Case(strcmp(q.selection.System,sys)));
    assert(isequal(chosen,actual(:)),'PMUAudit:Selection','QP selection differs from frozen protocol.');
    for k=1:45
        inst=d.instances{k};G=sparse(inst.G);m=inst.m;truth=inst.true_xi(:);
        ys=[inst.ytrue(:)+inst.noise(:),inst.y(:)];
        assert(norm(ys(:,2)-ys(:,1)-inst.bad(:),inf)<=1e-12*max(1,norm(ys(:,2),inf)));
        assert(norm(G*truth-inst.ytrue(:))<=1e-10*max(1,norm(inst.ytrue)));
        % The six independent recoveries share one fresh sparse-LS call.
        indices=zeros(2,3);rhs=zeros(m,6);col=0;
        for sc=1:2
            if sc==1,scenario='GaussianOnly';else,scenario='SparseBadData';end
            for a=1:3
                ix=find(strcmp({w.rows.System},sys)&[w.rows.Case]==k ...
                    &strcmp({w.rows.Scenario},scenario)&strcmp({w.rows.Method},names{a}));
                assert(isscalar(ix)&&~checked(ix),'PMUAudit:Duplicate','Missing/duplicate study row.');
                indices(sc,a)=ix;col=col+1;r=w.endpoints{ix};
                if a==1,e=zeros(m,1);else,e=split_error(r,m);end
                rhs(:,col)=ys(:,sc)-e;
            end
        end
        good=all(isfinite(rhs),1);states=nan(size(G,2),6);
        if any(good),states(:,good)=G\rhs(:,good);end
        col=0;
        for sc=1:2
            y=ys(:,sc);
            for a=1:3
                col=col+1;ix=indices(sc,a);row=w.rows(ix);r=w.endpoints{ix};xi=states(:,col);
                assert(row.Currents==inst.Currents&&row.BadDataIfContaminated==nnz(inst.bad) ...
                    &&row.RandomSeed==inst.RandomSeed&&row.Measurements==m);
                if a==1
                    fit=y-G*xi;normal=norm(G'*fit,inf)/max(1,norm(G,1)* ...
                        (norm(y,inf)+norm(G,inf)*norm(xi,inf)));
                    assert(all(isfinite(xi))&&normal<=1e-12&&isnan(row.RobustSolved), ...
                        'PMUAudit:WLS','Invalid LS solution or WLS labelled with robust success.');
                    re=norm(xi-truth)/norm(truth);
                    close_value(re,r.relerr,1e-10,'WLS endpoint state error');
                    close_value(re,row.RelativeStateError,1e-10,'WLS row state error');
                    close_value(normal,row.LSNormalResidualScaled,1e-12,'WLS normal residual');
                    sd=norm(xi-r.xi);assert(sd<=1e-10*max(1,norm(xi)));
                    stats.WLSChecks=stats.WLSChecks+1;
                    stats.MaxWLSNormalDataScaled=max(stats.MaxWLSNormalDataScaled,normal);
                    stats.MaxFreshStateDifference=max(stats.MaxFreshStateDifference,sd);
                    fresh_re(ix)=re;fresh_success(ix)=true;
                else
                    z=robust_metrics(r,G,y,truth,c,xi);stats=update_stats(stats,z,r);
                    if z.Finite
                        close_value(row.OriginalResidual,z.Residual,1e-9,'WLS-study row residual');
                        close_value(row.Objective,z.Objective,1e-10,'WLS-study row objective');
                        close_value(row.RelativeStateError,z.RE,1e-10,'WLS-study row state error');
                    else,assert(~r.success&&~logical(row.RobustSolved));end
                    assert(logical(row.RobustSolved)==logical(r.success));
                    if a==2
                        check_nf(r);stats.NativeNFChecks=stats.NativeNFChecks+1;
                    else,assert(r.monotonicity_violations==0);end
                    fresh_re(ix)=z.RE;fresh_success(ix)=logical(r.success);
                    if formal_available&&sc==2&&~strcmp(sys,'ieee300')
                        fi=find(strcmp({f.rows.System},sys)&[f.rows.Case]==k ...
                            &strcmp({f.rows.Method},names{a})&[f.rows.Repeat]==1);
                        assert(isscalar(fi),'PMUAudit:FormalCoverage','Missing formal comparison record.');
                        old=f.endpoints{fi};
                        same_backend=strcmp(system_backend(f.protocol.Config,sys), ...
                            system_backend(w.protocol.Config,sys));
                        if same_backend
                            assert(logical(old.success)==logical(r.success));
                            if a==2
                                assert(r.ActualItr==old.ActualItr&&r.ActualNF==old.ActualNF);
                            else,assert(r.Itr==old.Itr);end
                            delta=norm(r.s-old.s)/max(1,norm(old.s));
                            assert(isfinite(delta)&&delta<=1e-10,'PMUAudit:Reproduction','Formal endpoint changed.');
                            stats.MaxFormalEndpointRelativeDifference=max(stats.MaxFormalEndpointRelativeDifference,delta);
                            stats.FormalExactReproductionChecks=stats.FormalExactReproductionChecks+1;
                        else
                            % Equivalent projectors can have different floating-point
                            % trajectories. Validate the reference endpoint by fresh LS.
                            robust_metrics(old,G,y,truth,c,[]);
                            stats.FormalCrossBackendChecks=stats.FormalCrossBackendChecks+1;
                        end
                        stats.FormalReproductionChecks=stats.FormalReproductionChecks+1;
                    end
                end
                close_value(row.RelativeStateErrorPercent,100*fresh_re(ix),1e-10,'Percent conversion');
                assert(logical(row.ValidState)==all(isfinite(r.xi)));
                checked(ix)=true;
            end
        end
        pi=find(strcmp({w.paired.System},sys)&[w.paired.Case]==k);assert(isscalar(pi));
        pr=w.paired(pi);fields={'WLSCleanRE','TSCleanRE','IRLSCleanRE';'WLSBadRE','TSBadRE','IRLSBadRE'};
        for sc=1:2
            for a=1:3,close_value(pr.(fields{sc,a}),fresh_re(indices(sc,a)),1e-10,'Paired state error');end
        end
        assert(logical(pr.CommonRobustSuccess)==all(fresh_success(indices(:,2:3)),'all'));
        stats.PairChecks=stats.PairChecks+1;
        qi=find(strcmp({q.rows.System},sys)&[q.rows.Case]==k);
        if ~isempty(qi)
            assert(isscalar(qi)&&~qp_checked(qi));
            sel=q.selection(qi,:);assert(sel.Case==k&&strcmp(sel.System{1},sys) ...
                &&sel.Measurements==m&&sel.BadData==nnz(inst.bad)&&sel.RandomSeed==inst.RandomSeed);
            a=q.endpoints{qi,1};b=q.endpoints{qi,2};
            za=robust_metrics(a,G,ys(:,2),truth,c,[]);zb=robust_metrics(b,G,ys(:,2),truth,c,[]);
            stats=update_stats(stats,za,a);stats=update_stats(stats,zb,b);
            check_nf(a);stats.NativeNFChecks=stats.NativeNFChecks+1;
            row=q.rows(qi);assert(row.TSVerified==a.success&&row.QPVerified==b.success);
            close_value(row.TSOriginalResidual,za.Residual,1e-9,'QP-study TS residual');
            close_value(row.QPOriginalResidual,zb.Residual,1e-9,'QP-study QP residual');
            close_value(row.TSObjective,za.Objective,1e-10,'QP-study TS objective');
            close_value(row.QPObjective,zb.Objective,1e-10,'QP-study QP objective');
            close_value(row.AbsoluteObjectiveDifference,abs(za.Objective-zb.Objective),1e-10,'Objective difference');
            close_value(row.RelativeObjectiveDifference,abs(za.Objective-zb.Objective) ...
                /max(abs(zb.Objective),realmin),1e-10,'Relative objective difference');
            state_distance=NaN;
            if za.Finite&&zb.Finite
                state_distance=norm(za.State-zb.State)/max(norm(zb.State),realmin);
            end
            close_value(row.RelativeStateVectorDifference,state_distance,1e-10,'Relative state-vector difference');
            stats.QPChecks=stats.QPChecks+1;
            if b.success
                assert(b.exitflag>0&&b.NumericKKTVerified&&b.NumericDualityVerified);
                [kk,gap]=qp_optimality(b,G,ys(:,2),c,zb.Objective);
                assert(kk<=c.target&&abs(gap)<=c.target,'PMUAudit:QPOptimality','QP KKT/duality check failed.');
                stats.QPVerified=stats.QPVerified+1;stats.QPKKTChecks=stats.QPKKTChecks+1;
                stats.QPDualityChecks=stats.QPDualityChecks+1;
                stats.MaxFreshQPKKTScaled=max(stats.MaxFreshQPKKTScaled,kk);
                stats.MaxFreshQPScaledDualityGap=max(stats.MaxFreshQPScaledDualityGap,abs(gap));
            else,stats.RetainedQPFailures=stats.RetainedQPFailures+1;end
            qp_checked(qi)=true;
        end
    end
end
assert(all(checked)&&all(qp_checked)&&stats.WLSChecks==270&&stats.PairChecks==135 ...
    &&stats.RobustChecks==560&&stats.NativeNFChecks==280&&stats.FormalReproductionChecks==180*formal_available);
check_summaries(w.summary,w.rows,w.paired);
assert(height(w.strata)==162,'PMUAudit:Strata','Expected 27 strata times six method/condition rows.');
for j=1:height(w.strata)
    st=w.strata(j,:);mask=strcmp({w.rows.System},st.System{1}) ...
        &[w.rows.Currents]==st.Currents&[w.rows.BadDataIfContaminated]==st.BadData;
    pp=w.paired(strcmp({w.paired.System},st.System{1}) ...
        &[w.paired.Currents]==st.Currents&[w.paired.BadData]==st.BadData);
    assert(st.Cases==5&&numel(pp)==5);check_summaries(st,w.rows(mask),pp);
end
report=struct('Created',char(datetime('now','TimeZone','Asia/Shanghai')),'Pass',true, ...
    'MATLAB',version,'ParameterProfile',c.profile,'InputCases',135, ...
    'WLSStudyEndpoints',810,'WLSStudyRobustEndpoints',540,'QPStudyEndpoints',20, ...
    'OriginalTolerance',c.target,'FeasibilityTolerance',c.feasibility_target, ...
    'SourceSnapshotChecks',source_checks,'FormalReferenceAvailable',formal_available,'Statistics',stats,'RuntimeRankingReported',false, ...
    'Scope','Supplemental paired WLS accuracy and ten-case QP reference; not full three-grid repeated CGPM timing.');
save(fullfile(out,'supplement_audit.mat'),'report');
fid=fopen(fullfile(out,'supplement_audit.json'),'w','n','UTF-8');assert(fid>=0);
filecleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));disp(report);
end

function e=split_error(r,m)
assert(isfield(r,'s')&&numel(r.s)==2*m,'PMUAudit:Split','Wrong endpoint dimension.');
s=r.s(:);e=s(1:m)-s(m+1:end);
end
function z=robust_metrics(r,G,y,truth,c,xi)
m=numel(y);
z=struct('Finite',false,'Residual',Inf,'Feasibility',Inf,'Objective',NaN,'RE',NaN, ...
    'ResidualDifference',0,'ObjectiveDifference',0,'REDifference',0,'StateDifference',0,'State',[]);
if ~isfield(r,'s')||numel(r.s)~=2*m,assert(~r.success);return;end
e=split_error(r,m);s=r.s(:);
if any(~isfinite(s)),assert(~r.success);return;end
if isempty(xi),xi=G\(y-e);end
fit=y-e-G*xi;fv=[min(s(1:m),c.lambda-fit);min(s(m+1:end),c.lambda+fit)];
z.Residual=norm(fv);z.Feasibility=norm(min(s,0));
z.Objective=.5*dot(fit,fit)+c.lambda*norm(e,1);z.RE=norm(xi-truth)/norm(truth);
z.Finite=all(isfinite([xi;fv;z.Objective;z.RE]));z.State=xi;
if ~z.Finite,assert(~r.success);return;end
close_value(z.Residual,r.residual,1e-9,'Native residual');
close_value(z.Objective,r.obj,1e-10,'Native objective');
close_value(z.RE,r.relerr,1e-10,'Native state error');
close_value(z.Feasibility,r.feasibility,1e-12,'Native feasibility');
assert(norm(e-r.e)<=1e-12*max(1,norm(e)));
z.StateDifference=norm(xi-r.xi);assert(z.StateDifference<=1e-10*max(1,norm(xi)));
z.ResidualDifference=abs(z.Residual-r.residual);z.ObjectiveDifference=abs(z.Objective-r.obj);
z.REDifference=abs(z.RE-r.relerr);
if r.success
    assert(z.Residual<=c.target&&z.Feasibility<=c.feasibility_target, ...
        'PMUAudit:FalseSuccess','Successful endpoint fails strict original tolerance.');
end
end
function stats=update_stats(stats,z,r)
stats.RobustChecks=stats.RobustChecks+1;
if r.success
    stats.RobustSuccesses=stats.RobustSuccesses+1;
    stats.MaxFreshSuccessfulResidual=max(stats.MaxFreshSuccessfulResidual,z.Residual);
else,stats.RetainedRobustFailures=stats.RetainedRobustFailures+1;end
if ~z.Finite,stats.NonfiniteFailedEndpoints=stats.NonfiniteFailedEndpoints+1;return;end
stats.MaxFreshResidualDifference=max(stats.MaxFreshResidualDifference,z.ResidualDifference);
stats.MaxFreshObjectiveDifference=max(stats.MaxFreshObjectiveDifference,z.ObjectiveDifference);
stats.MaxFreshRelativeErrorDifference=max(stats.MaxFreshRelativeErrorDifference,z.REDifference);
stats.MaxFreshStateDifference=max(stats.MaxFreshStateDifference,z.StateDifference);
end
function check_nf(r)
assert(r.count_ok&&r.report_ok&&r.NFCounted==r.ActualNF&&r.NFCounted==r.info.actual_NF, ...
    'PMUAudit:NF','Native and independently counted NF disagree.');
end
function [kk,scaled_gap]=qp_optimality(r,G,y,c,obj)
m=numel(y);n=size(G,2);x=r.QPSolution(:);assert(numel(x)==n+3*m);
xi=x(1:n);e=x(n+(1:m));v=x(n+m+(1:m));t=x(n+2*m+(1:m));
nu=r.QPMultipliers.eqlin(:);mu=r.QPMultipliers.ineqlin(:);assert(numel(nu)==m&&numel(mu)==2*m);
mp=mu(1:m);mn=mu(m+1:end);eq=G*xi+e+v-y;ineq=[e-t;-e-t];
grad=[zeros(n+m,1);v;c.lambda*ones(m,1)];
di=[zeros(n,1);mp-mn;zeros(m,1);-mp-mn];de=[G'*nu;nu;nu;zeros(m,1)];
native_obj=.5*dot(v,v)+c.lambda*sum(t);
kk=max([norm(eq,inf)/max([1,norm(y,inf),norm(G*xi+e+v,inf)]), ...
    max([0;ineq])/max(1,norm(ineq,inf)),max([0;-mu])/max(1,norm(mu,inf)), ...
    norm(grad+di+de,inf)/max([1,norm(grad,inf),norm(di,inf),norm(de,inf)]), ...
    norm(mu.*(-ineq),inf)/max(1,abs(native_obj))]);
close_value(kk,r.KKT.ResidualScaled,1e-10,'QP KKT residual');
dual=r.Duality.Vector(:);assert(numel(dual)==m);
assert(norm(G'*dual,inf)/max(1,norm(G,1)*norm(dual,inf))<=c.target ...
    &&max(0,norm(dual,inf)-c.lambda)<=c.target,'PMUAudit:QPDual','Invalid numerical dual point.');
dual_obj=dot(y,dual)-.5*dot(dual,dual);gap=obj-dual_obj;
scaled_gap=gap/max([1,abs(obj),abs(dual_obj)]);
close_value(gap,r.Duality.Gap,1e-10,'QP duality gap');
end
function check_summaries(S,rows,paired)
for j=1:height(S)
    z=S(j,:);rr=rows(strcmp({rows.System},z.System{1}) ...
        &strcmp({rows.Scenario},z.Scenario{1})&strcmp({rows.Method},z.Method{1}));
    pp=paired(strcmp({paired.System},z.System{1}));ids=[pp([pp.CommonRobustSuccess]).Case];
    v=[rr.RelativeStateErrorPercent];v=v(isfinite(v));qs=prctile(v,[25 50 75]);
    cv=[rr(ismember([rr.Case],ids)).RelativeStateErrorPercent];
    assert(z.Cases==numel(rr)&&z.ValidStates==sum([rr.ValidState])&&z.CommonBothScenariosCases==numel(ids));
    if strcmp(z.Method{1},'WLS'),assert(isnan(z.RobustSolved));
    else,assert(z.RobustSolved==sum([rr.RobustSolved]));end
    close_value(z.MeanREPercentAll,mean(v),1e-10,'Summary mean');
    close_value(z.MedianREPercentAll,qs(2),1e-10,'Summary median');
    close_value(z.Q1REPercentAll,qs(1),1e-10,'Summary Q1');
    close_value(z.Q3REPercentAll,qs(3),1e-10,'Summary Q3');
    close_value(z.MeanREPercentCommon,mean(cv),1e-10,'Common-pair summary mean');
end
end
function check_config(a,b)
fields={'lambda','target','feasibility_target','budget','backend','initial','profile'};
for j=1:numel(fields),assert(isequal(a.(fields{j}),b.(fields{j})),'PMUAudit:Config','Configuration changed.');end
for sys=b.systems
    assert(strcmp(system_backend(a,sys{1}),system_backend(b,sys{1})), ...
        'PMUAudit:Backend','Recorded per-system projector differs.');
end
p=a.params{1};q=b.params{1};assert(isequal(fieldnames(p),fieldnames(q)));
for f=fieldnames(p).'
    v=p.(f{1});u=q.(f{1});
    if isa(v,'function_handle'),assert(isa(u,'function_handle')&&strcmp(func2str(v),func2str(u)));
    else,assert(isequaln(v,u),'PMUAudit:TSParameter','TS parameter changed.');end
end
end
function backend=system_backend(c,sys)
backend=c.backend;
if isfield(c,'backend_by_system')&&isfield(c.backend_by_system,sys)
    backend=c.backend_by_system.(sys);
end
end
function n=check_sources(snapshot,names,here)
n=0;
for j=1:numel(names)
    ix=find(strcmp({snapshot.Name},names{j}));assert(isscalar(ix),'PMUAudit:Source','Source absent from snapshot.');
    assert(strcmp(snapshot(ix).Text,fileread(fullfile(here,names{j}))), ...
        'PMUAudit:SourceChanged','Core source changed since the saved study: %s',names{j});
    n=n+1;
end
end
function close_value(a,b,tol,label)
if isnan(a)&&isnan(b),return;end
if isinf(a)&&isequal(a,b),return;end
assert(isscalar(a)&&isscalar(b)&&isfinite(a)&&isfinite(b) ...
    &&abs(a-b)<=tol*max([1,abs(a),abs(b)]),'PMUAudit:Mismatch','Mismatch: %s.',label);
end
function h=sha256_file(p)
fid=fopen(p,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),b=fread(fid,1048576,'*uint8');if ~isempty(b),md.update(typecast(b,'int8'));end;end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function restore_session(p,t)
path(p);maxNumCompThreads(t);
end
