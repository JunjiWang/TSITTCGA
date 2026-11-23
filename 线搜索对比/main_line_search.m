function main_line_search(options)
% Controlled acceptance-criterion comparison, then warm-start ablation.
% main_line_search
% main_line_search(struct('verify_only',true))
% main_line_search(struct('figs_only',true))
if nargin==0, options=struct; end
if isfield(options,'mode') && strcmp(options.mode,'step_probe')
    step_probe_study(rmfield(options,'mode')); return
end
if isfield(options,'mode') && strcmp(options.mode,'initial_steps')
    initial_step_study(rmfield(options,'mode')); return
end
if isfield(options,'fair') && options.fair
    fair_study(rmfield(options,'fair')); return
end
if ~isfield(options,'legacy') || ~options.legacy
    abc_study(options); return
end
options=rmfield(options,'legacy');
o=struct('verify_only',false,'figs_only',false,'plots',true,'resume',false, ...
    'groups',{{'monotone','pseudo','signal'}},'phases',{{'criteria','warm_start','criteria_sigma010','criteria_sigma050'}});
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'Unknown option: %s',key{1}); o.(key{1})=options.(key{1});
end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
outdir=fullfile(here,'results'); figdir=fullfile(here,'figures');
savedpath=path; savedrng=rng; cleanup=onCleanup(@()restore_session(savedpath,savedrng)); %#ok<NASGU>
addpath(here,'-begin');
if o.figs_only
    r=load(fullfile(outdir,'experiments.mat'),'experiments');
    write_tables(r.experiments,outdir); plot_all(r.experiments,figdir); return;
end
verify_core(root,here);
if o.verify_only, fprintf('VERIFICATION_PASSED\n'); return; end
if ~isfolder(outdir), mkdir(outdir); end
if o.plots && ~isfolder(figdir), mkdir(figdir); end
experiments=struct;
protocol=struct('criteria_sigma',.01,'criteria_initial',1,'criteria_shrink',.5, ...
    'criterion_names',{{'Proposed','Yin2021','Liu','Jiang'}}, ...
    'warm_start_flags',[false false true true],'increase_flags',[false true false true], ...
    'function_tolerance',1e-6,'feasibility_tolerance',1e-6,'iteration_limit',2000, ...
    'source_project',root,'uses_formal_initial_points',true,'parameter_tuning',false);
source_files={'main_line_search.m','TSITTCGA_LS.m'};
source_snapshot=cellfun(@(q)fileread(fullfile(here,q)),source_files,'UniformOutput',false);
matlab_version=version; matlab_threads=maxNumCompThreads; timer=tic;
if o.resume && isfile(fullfile(outdir,'experiments.mat'))
    previous=load(fullfile(outdir,'experiments.mat'),'experiments','protocol','source_snapshot');
    assert(isequaln(previous.protocol,protocol),'Saved protocol differs from current protocol.');
    assert(strcmp(previous.source_snapshot{2},source_snapshot{2}),'Solver changed since saved computations.');
    experiments=previous.experiments;
    for fld=fieldnames(experiments)'
        if ~isfield(experiments.(fld{1}),'execution_sources')
            experiments.(fld{1}).execution_sources=previous.source_snapshot;
        end
    end
end
for ph=1:numel(o.phases)
    phase=o.phases{ph}; assert(any(strcmp(phase,{'criteria','warm_start','criteria_sigma010','criteria_sigma050'})));
    for g=1:numel(o.groups)
        group=o.groups{g};
        fieldkey=[phase,'_',group];
        if isfield(experiments,fieldkey)
            cached=experiments.(fieldkey); src=source_dir(root,group);
            assert(strcmp(cached.input_meta.original_main,fileread(fullfile(src,'main.m'))));
            assert(strcmp(cached.input_meta.original_solver,fileread(fullfile(src,'TSITTCGA.m'))));
            pf='problems.m'; if strcmp(group,'signal'), pf='problem.m'; end
            assert(strcmp(cached.input_meta.problem_source,fileread(fullfile(src,pf))));
            fprintf('REUSE saved %s / %s\n',phase,group);
            if o.plots && ~startsWith(phase,'criteria_sigma'), plot_result(cached,figdir); end
            print_summary(cached); continue;
        end
        fprintf('\nSTART %s / %s\n',phase,group);
        result=run_group(root,group,phase); result.execution_sources=source_snapshot; experiments.([phase,'_',group])=result;
        elapsed=toc(timer);
        save(fullfile(outdir,'experiments.mat'),'experiments','protocol','o', ...
            'source_files','source_snapshot','matlab_version','matlab_threads','elapsed','-v7.3');
        write_tables(experiments,outdir);
        if o.plots && ~startsWith(phase,'criteria_sigma'), plot_result(result,figdir); end
        print_summary(result);
        fprintf('SAVED %s / %s. Elapsed %.1f s\n',phase,group,elapsed);
    end
end
fprintf('\nLINE_SEARCH_EXPERIMENTS_DONE %.1f s\n',toc(timer));
end

function verify_core(root,here,groups)
if nargin<3, groups={'monotone','pseudo','signal'}; end
q=struct('Itr_max',1,'sigma',.5,'mu',.8,'kappa',.5,'a',.01,'rho_ls',.8, ...
    'varsigma',1,'rho',1,'eta',.1,'tbar',.2,'c',.05,'r',.05, ...
    'epsilon_seq',@(j)1/(j+1)^2,'quiet',true);
pr=struct('name','CriterionCheck','f',@identity_problem);
cfg=make_configs(q,'criteria'); expected=[.4 .8 .8 .8];
for j=1:4
    c=cfg{j}; c.sigma=.5; c.initial=.8;
    [~,~,~,~,info]=TSITTCGA_LS(pr,1,q,c);
    assert(abs(info.accepted_steps(1)-expected(j))<1e-14,'Unexpected first accepted step.');
    assert(info.actual_NF==1+info.inertia_nf+info.ls_nf+info.projection_nf);
end
bad=struct('name','InfeasibleTrialRoot','f',@infeasible_root_problem);
for j=1:4
    [~,~,~,~,info]=TSITTCGA_LS(bad,1,q,cfg{j});
    assert(~info.success && info.backtracks>=1 && info.accepted_steps(1)<1);
end
for group=groups
    group=group{1}; src=source_dir(root,group); saved=path; cl=onCleanup(@()path(saved)); %#ok<NASGU>
    addpath(src,'-begin'); addpath(here,'-begin');
    [cases,p]=load_cases(root,group); f=str2func('TSITTCGA');
    assert(strcmp(functions(f).file,fullfile(src,'TSITTCGA.m')),'Wrong original solver.');
    configs=make_configs(p,'warm_start');
    for i=unique([1,numel(cases)])
        short=p; short.Itr_max=8;
        [~,~,~,~,a]=f(cases(i).PR,cases(i).x0,short);
        [~,~,~,~,b]=TSITTCGA_LS(cases(i).PR,cases(i).x0,short,configs{4});
        assert(isequaln(a.x,b.x) && a.actual_NF==b.actual_NF ...
            && a.loop_iterations==b.loop_iterations && a.success==b.success, ...
            'Native-core mismatch for %s.',group);
    end
    clear cl;
end
fprintf('Criterion equations, root safeguards, NF accounting and native-core parity checked.\n');
end

function out=identity_problem(x,mode)
out=x; %#ok<INUSD>
end
function out=infeasible_root_problem(x,mode)
if mode==1, out=x+1; else, out=max(x,0); end
end

function result=run_group(root,group,phase)
src=source_dir(root,group); saved=path; cl=onCleanup(@()path(saved)); %#ok<NASGU>
addpath(src,'-begin'); [cases,p,meta]=load_cases(root,group);
cfg=make_configs(p,phase); nc=numel(cases); nm=numel(cfg);
reps=1; if strcmp(group,'pseudo'), reps=5; elseif strcmp(group,'signal'), reps=3; end
M=nan(nc,nm); result=struct('phase',phase,'group',group,'params',p, ...
    'configurations',{cfg},'labels',{cellfun(@(q)q.label,cfg,'UniformOutput',false)}, ...
    'cases',rmfield(cases,{'PR','x0'}),'input_meta',meta,'timing_repeats',reps, ...
    'Itr',M,'NF',M,'LSNF',M,'Time',M,'NormF',M,'Feasibility',M,'MSE',M, ...
    'Backtracks',M,'IncreaseTests',M,'IncreaseAccepts',M,'InitialAccepts',M, ...
    'MedianStep',M,'Success',false(nc,nm),'Status',{cell(nc,nm)}, ...
    'RawTime',nan(nc,nm,reps),'RawNF',nan(nc,nm,reps));
short=p; short.Itr_max=3;
for j=1:nm, TSITTCGA_LS(cases(1).PR,cases(1).x0,short,cfg{j}); end
wall=tic;
for i=1:nc
    first=cell(1,nm);
    for r=1:reps
        for j=mod(i+r-2+(0:nm-1),nm)+1
            [~,nf,it,fn,d]=TSITTCGA_LS(cases(i).PR,cases(i).x0,p,cfg{j});
            fx=cases(i).PR.f(d.x,1); checked=norm(fx);
            feasible=norm(d.x-cases(i).PR.f(d.x,2));
            good=d.success && isreal(d.x) && all(isfinite(d.x)) && isreal(fx) ...
                && isfinite(checked) && checked<=1e-6 && isfinite(feasible) ...
                && feasible<=1e-6 && it<=p.Itr_max;
            assert(~d.success || good,'Invalid reported success.');
            assert(d.actual_NF==1+d.inertia_nf+d.ls_nf+d.projection_nf);
            signature={it,nf,fn,checked,feasible,d.success,d.status,d.x, ...
                d.ls_nf,d.backtracks,d.increase_tests,d.increase_accepts,d.accepted_steps};
            if r==1
                first{j}=signature; result.Success(i,j)=good; result.Status{i,j}=d.status;
                result.Itr(i,j)=it; result.NF(i,j)=nf; result.LSNF(i,j)=d.ls_nf;
                result.NormF(i,j)=checked; result.Feasibility(i,j)=feasible;
                result.Backtracks(i,j)=d.backtracks; result.IncreaseTests(i,j)=d.increase_tests;
                result.IncreaseAccepts(i,j)=d.increase_accepts; result.InitialAccepts(i,j)=d.initial_accepts;
                if ~isempty(d.accepted_steps), result.MedianStep(i,j)=median(d.accepted_steps); end
                if strcmp(group,'signal'), result.MSE(i,j)=d.mse(end); end
            else
                assert(isequaln(first{j},signature),'Nondeterministic costs/iterates.');
            end
            result.RawTime(i,j,r)=d.actual_time; result.RawNF(i,j,r)=d.actual_NF;
        end
    end
    for j=1:nm
        if result.Success(i,j), result.Time(i,j)=median(result.RawTime(i,j,:)); end
    end
    if ~strcmp(group,'monotone') || mod(i,5)==0
        fprintf('%s / %s: %d/%d cases, %.1f s\n',phase,group,i,nc,toc(wall));
    end
end
result.wall=toc(wall);
end

function cfg=make_configs(p,phase)
c=struct('criterion','Proposed','sigma',.01,'initial',1,'shrink',.5, ...
    'warm',false,'increase',false,'lambda',.001,'nu',.8,'liu_mu',.01, ...
    'j_a',2,'j_l',.5,'label','Proposed');
cfg=repmat({c},1,4);
if startsWith(phase,'criteria')
    keys={'Proposed','Yin2021','Liu','Jiang'}; labels={'Proposed (this paper)','Yin et al. [10]','Liu [15]','Jiang et al. [16]'};
    sig=.01;
    if strcmp(phase,'criteria_sigma010'), sig=.1; elseif strcmp(phase,'criteria_sigma050'), sig=.5; end
    for j=1:4
        cfg{j}.criterion=keys{j}; cfg{j}.label=labels{j}; cfg{j}.sigma=sig;
    end
else
    labels={'Fixed','Fixed + increase','Warm start','Warm start + increase'};
    for j=1:4
        cfg{j}.sigma=p.sigma; cfg{j}.initial=p.rho_ls; cfg{j}.shrink=p.kappa;
        cfg{j}.warm=j>=3; cfg{j}.increase=mod(j,2)==0; cfg{j}.label=labels{j};
    end
end
end

function [cases,p,meta]=load_cases(root,group)
src=source_dir(root,group); file=fullfile(src,'results','compare.mat');
cases=struct('Problem',{},'n',{},'Start',{},'PR',{},'x0',{});
if strcmp(group,'monotone')
    s=load(file,'wjj','inst','X0'); p=s.wjj;
    clear problems; pf=str2func('problems');
    assert(strcmp(functions(pf).file,fullfile(src,'problems.m')),'Wrong monotone problem.');
    for i=1:numel(s.inst)
        id=s.inst(i).id; pr=struct('name',sprintf('P%d',id),'n',s.inst(i).n,'f',@(x,mode)pf(x,mode,id));
        for r=1:size(s.X0,2)
            cases(end+1)=struct('Problem',pr.name,'n',pr.n,'Start',r,'PR',pr,'x0',s.X0{i,r}); %#ok<AGROW>
        end
    end
    meta=struct('original_result',file,'X0',{s.X0},'problem_ids',[s.inst.id],'dimensions',[s.inst.n]);
elseif strcmp(group,'pseudo')
    s=load(file,'wjj','X0','problem'); p=s.wjj;
    clear problems; pf=str2func('problems');
    assert(strcmp(functions(pf).file,fullfile(src,'problems.m')),'Wrong pseudo-monotone problem.');
    for i=1:numel(s.X0)
        id=s.problem(i); label='P11'; r=i-8;
        if id==9, label='P9R'; r=i; end
        pr=struct('name',label,'n',numel(s.X0{i}),'f',@(x,mode)pf(x,mode,id));
        cases(end+1)=struct('Problem',label,'n',pr.n,'Start',r,'PR',pr,'x0',s.X0{i}); %#ok<AGROW>
    end
    meta=struct('original_result',file,'X0',{s.X0},'problem_ids',s.problem);
elseif strcmp(group,'signal')
    s=load(file,'p1','f','x0','n','m','n_spikes','noise_sigma','lambda'); p=s.p1;
    rng(2016,'twister'); rand;
    n=s.n; m=s.m; f=zeros(n,1); q=randperm(n); f(q(1:s.n_spikes))=randn(s.n_spikes,1);
    A=randn(m,n)/sqrt(n); AT=A'; b=A*f+s.noise_sigma*randn(m,1);
    Atb=AT*b; lambda=.005*max(abs(Atb)); x0=[max(Atb,0);max(-Atb,0)];
    assert(isequal(f,s.f) && isequal(x0,s.x0) && isequal(lambda,s.lambda),'Signal-data mismatch.');
    clear problem; pf=str2func('problem');
    assert(strcmp(functions(pf).file,fullfile(src,'problem.m')),'Wrong signal problem.');
    pr=struct('name','SparseSignalRecovery','n',n,'signal_n',n,'true_x',f, ...
        'f',@(x,mode)pf(x,mode,A,AT,Atb,lambda,n));
    cases=struct('Problem','SparseSignalRecovery','n',n,'Start',1,'PR',pr,'x0',x0);
    meta=struct('original_result',file,'rng_seed',2016,'rng_initial_rand',true, ...
        'f',f,'x0',x0,'n',n,'m',m,'n_spikes',s.n_spikes,'noise_sigma',s.noise_sigma,'lambda',lambda);
else
    error('Unknown group.');
end
assert(p.Itr_max==2000,'Unexpected original iteration limit.');
meta.original_main=fileread(fullfile(src,'main.m'));
problem_file='problems.m'; if strcmp(group,'signal'), problem_file='problem.m'; end
meta.problem_source=fileread(fullfile(src,problem_file));
meta.original_solver=fileread(fullfile(src,'TSITTCGA.m'));
end

function src=source_dir(root,group)
switch group
    case 'monotone', src=fullfile(root,'单调方程组');
    case 'pseudo', src=fullfile(root,'伪单调方程组');
    case 'signal', src=fullfile(root,'信号恢复');
    otherwise, error('Unknown group.');
end
end

function write_tables(experiments,outdir)
rows=cell(0,19); summary=cell(0,13); comparisons=cell(0,13);
for field=fieldnames(experiments)'
    s=experiments.(field{1}); nc=size(s.NF,1);
    for j=1:4
        for i=1:nc
            rows(end+1,:)={s.phase,s.group,s.cases(i).Problem,s.cases(i).n,s.cases(i).Start, ...
                s.labels{j},s.Itr(i,j),s.NF(i,j),s.LSNF(i,j),s.Time(i,j),s.NormF(i,j), ...
                s.Feasibility(i,j),s.Backtracks(i,j),s.IncreaseTests(i,j),s.IncreaseAccepts(i,j), ...
                s.MedianStep(i,j),s.MSE(i,j),s.Success(i,j),s.Status{i,j}}; %#ok<AGROW>
        end
        % Full-set means, not separate deletion of difficult cases per variant.
        summary(end+1,:)={s.phase,s.group,s.labels{j},nc,mean(s.Success(:,j)), ...
            mean(s.Itr(:,j)),mean(s.NF(:,j)),complete_mean(s.LSNF(:,j),s.Success(:,j)),mean(s.Time(:,j)), ...
            complete_mean(s.Backtracks(:,j),s.Success(:,j)),complete_mean(s.IncreaseTests(:,j),s.Success(:,j)),complete_mean(s.IncreaseAccepts(:,j),s.Success(:,j)), ...
            max(s.NormF(:,j))}; %#ok<AGROW>
    end
    pairs=[1 2;1 3;1 4]; if strcmp(s.phase,'warm_start'), pairs=[3 1;4 2;2 1;4 3;4 1]; end
    for k=1:size(pairs,1)
        a=pairs(k,1); b=pairs(k,2); mask=s.Success(:,a)&s.Success(:,b);
        na=s.NF(mask,a); nb=s.NF(mask,b);
        comparisons(end+1,:)={s.phase,s.group,[s.labels{a},' / ',s.labels{b}],sum(mask), ...
            paired_ratio(s.NF(mask,a),s.NF(mask,b)),paired_ratio(s.LSNF(mask,a),s.LSNF(mask,b)), ...
            paired_ratio(s.Itr(mask,a),s.Itr(mask,b)),paired_ratio(s.Time(mask,a),s.Time(mask,b)), ...
            sum(na<nb),sum(na==nb),sum(na>nb),mean(s.Success(:,a)),mean(s.Success(:,b))}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'Phase','Group','Problem','n','Start','Variant','Itr','NF','LSNF', ...
    'Time','NormF','Feasibility','Backtracks','IncreaseTests','IncreaseAccepts','MedianStep', ...
    'MSE','Success','Status'});
writetable(T,fullfile(outdir,'all_runs.csv'),'Encoding','UTF-8');
T=cell2table(summary,'VariableNames',{'Phase','Group','Variant','Instances','SuccessRate', ...
    'MeanItr','MeanNF','MeanLSNF','MeanTime','MeanBacktracks','MeanIncreaseTests','MeanIncreaseAccepts','MaxNormF'});
writetable(T,fullfile(outdir,'summary.csv'),'Encoding','UTF-8');
T=cell2table(comparisons,'VariableNames',{'Phase','Group','Comparison','CommonInstances','NF_Ratio', ...
    'LSNF_Ratio','Itr_Ratio','Time_Ratio','NF_Wins','NF_Ties','NF_Losses','SuccessRateA','SuccessRateB'});
writetable(T,fullfile(outdir,'paired_comparisons.csv'),'Encoding','UTF-8');
end

function v=complete_mean(x,ok)
if all(ok), v=mean(x); else, v=NaN; end
end
function ratio=paired_ratio(a,b)
mask=isfinite(a)&isfinite(b)&a>0&b>0;
if any(mask), ratio=exp(mean(log(a(mask)./b(mask)))); else, ratio=NaN; end
end
function print_summary(s)
fprintf('Variant                         Mean Itr     Mean NF    LS NF      Time (s)\n');
for j=1:4
    fprintf('%-30s %10.2f %11.2f %9.2f %13.6f\n',s.labels{j},mean(s.Itr(:,j)), ...
        mean(s.NF(:,j)),complete_mean(s.LSNF(:,j),s.Success(:,j)),mean(s.Time(:,j)));
end
end

function plot_all(experiments,figdir)
if ~isfolder(figdir), mkdir(figdir); end
for field=fieldnames(experiments)'
    s=experiments.(field{1});
    if ~startsWith(s.phase,'criteria_sigma'), plot_result(s,figdir); end
end
end
function plot_result(s,figdir)
if ~isfolder(figdir), mkdir(figdir); end
blue=[0 114 178]/255; orange=[213 94 0]/255; green=[0 158 115]/255;
pink=[204 121 167]/255; gray=[90 90 90]/255;
colors=[blue;orange;green;pink]; if strcmp(s.phase,'warm_start'), colors=[gray;orange;green;blue]; end
styles={'-','--','-.',':'}; markers={'o','s','^','d'};
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 7.5]);
cl=onCleanup(@()close(f)); %#ok<NASGU>
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
fields={'NF','Itr'}; titles={'Function evaluations','Iterations'};
for m=1:2
    ax=nexttile; hold(ax,'on');
    identical=all(s.Success,'all') && all(s.NF==s.NF(:,1),'all') && all(s.Itr==s.Itr(:,1),'all');
    if strcmp(s.group,'signal') || identical
        values=mean(s.(fields{m}),1);
        vf=values(isfinite(values)); ymax=1; if ~isempty(vf), ymax=max(1,max(vf)); end
        for j=1:4
            bar(ax,j,values(j),.64,'FaceColor',colors(j,:),'EdgeColor','none');
            if isfinite(values(j))
                text(ax,j,values(j)+.025*ymax,sprintf('%g',round(values(j),1)), ...
                    'HorizontalAlignment','center','FontSize',10,'FontName','Times New Roman');
            end
        end
        set(ax,'XTick',1:4,'XTickLabel',s.labels,'XTickLabelRotation',18);
        ylim(ax,[0,ymax*1.25]);
        yl=titles{m}; if size(s.NF,1)>1, yl=['Mean ',lower(yl)]; end
        ylabel(ax,yl);
    else
        cost=s.(fields{m}); cost(~s.Success)=Inf; best=min(cost,[],2); ratios=cost./best;
        ratios(~isfinite(ratios))=Inf; finite=ratios(isfinite(ratios));
        upper=2; if ~isempty(finite), upper=max(2,max(finite)*1.03); end
        grid=unique([1;finite(:);upper]);
        for j=1:4
            y=arrayfun(@(q)mean(ratios(:,j)<=q*(1+1e-12)),grid);
            stairs(ax,log2(grid),y,'LineWidth',1.5,'Color',colors(j,:),'LineStyle',styles{j});
            ix=unique(round(linspace(1,numel(grid),min(5,numel(grid)))));
            plot(ax,log2(grid(ix)),y(ix),markers{j},'Color',colors(j,:),'MarkerSize',4, ...
                'LineStyle','none','HandleVisibility','off');
        end
        alphaUpper=max(1,ceil(log2(upper)));
        set(ax,'XScale','linear','XTick',unique([0:max(1,ceil(alphaUpper/6)):alphaUpper,alphaUpper]));
        xlim(ax,[0,alphaUpper]); ylim(ax,[0,1.025]);
        xlabel(ax,'\alpha = log_2(\tau)'); ylabel(ax,'\rho_s(\alpha)');
        lbl=s.labels;
lbl=strrep(lbl,'Proposed','Proposed (this paper)');
lbl=strrep(lbl,'Yin (2021)','Yin [2021]');
lbl=strrep(lbl,'Liu (2026)','Liu [2026]');
lbl=strrep(lbl,'Jiang et al. (2026)','Jiang [2026]');
legend(ax,lbl,'Location','southeast','Box','off','FontSize',7.2);
    end
    title(ax,titles{m},'FontWeight','normal');
    set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',.8, ...
        'TickDir','out','Box','on','YGrid','on','GridAlpha',.13);
end
stem=[s.phase,'_',s.group];
exportgraphics(f,fullfile(figdir,[stem,'.pdf']),'ContentType','vector');
exportgraphics(f,fullfile(figdir,[stem,'.png']),'Resolution',250);
end
function restore_session(savedpath,savedrng)
path(savedpath); rng(savedrng);
end

function abc_study(options)
% A: criterion only; B: fixed reduction versus warm start plus increase; C: complete algorithms.
o=struct('phases',{{'A','B','C'}},'verify_only',false,'figs_only',false, ...
    'resume',false,'plots',true,'threads',1,'output_dir','');
for k=fieldnames(options)'
    assert(isfield(o,k{1}),'Unknown ABC option %s',k{1}); o.(k{1})=options.(k{1});
end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
outbase=here;
if ~isempty(o.output_dir), outbase=o.output_dir; end
outdir=fullfile(outbase,'results'); figdir=fullfile(outbase,'figures');
if ~isfolder(outdir), mkdir(outdir); end
if ~isfolder(figdir), mkdir(figdir); end
outfile=fullfile(outdir,'experiments_ABC.mat');
if ~isempty(o.output_dir), outfile=fullfile(outdir,'compare.mat'); end
validateattributes(o.threads,{'numeric'},{'scalar','integer','positive','finite'});
assert(iscell(o.phases) && ~isempty(o.phases) && numel(unique(o.phases))==numel(o.phases), ...
    'phases must be a nonempty cell array without duplicate groups.');
if o.figs_only
    d=load(outfile,'experiments'); abc_plots(d.experiments,figdir); return
end
oldpath=path; oldrng=rng; oldthreads=maxNumCompThreads;
cleanup=onCleanup(@()abc_restore(oldpath,oldrng,oldthreads)); %#ok<NASGU>
maxNumCompThreads(o.threads);
src=source_dir(root,'monotone'); addpath(src,'-begin'); addpath(here,'-begin');
[cases,p,meta]=load_cases(root,'monotone');
native=load(fullfile(src,'results','compare.mat'),'names','params','source_files','source_snapshot');
assert(abc_same_value(p,abc_main_parameters(meta.original_main,p.Itr_max)), ...
    'Frozen TSITTCGA parameters do not match the current main.m.');
ix=find(strcmp(native.names,'TSITTCGA'));
assert(isscalar(ix) && abc_same_value(native.params{ix},p), ...
    'Formal TSITTCGA parameter metadata are inconsistent.');
assert(numel(cases)==200 && numel(unique([cases.n]))>=5,'Unexpected eight-problem design.');
assert(abs(p.sigma-0.80)<1e-9 && abs(p.rho_ls-0.80)<1e-9 && abs(p.kappa-0.50)<1e-9,'A control constants must follow frozen parameters.');
protocol=struct('date','2026-10-15','instances',200,'problem_count',8, ...
    'same_inputs',true,'max_iterations',2000,'epsilon',1e-6,'epsilon_C',1e-6, ...
    'A_sigma',p.sigma,'A_initial',p.rho_ls,'A_shrink',p.kappa, ...
    'A_warm',false,'A_increase',false,'B_increase_flags',[false true],'B_warm_flags',[false true], ...
    'A_criteria',{{'Proposed','Yin2021','Liu','Jiang'}}, ...
    'C_names',{{'TSITTCGA','AT2DFP','ADCGPM','TSIM-DFPM-PRP','HTTCGP'}}, ...
    'threads',o.threads,'uses_formal_starts',true,'parameter_tuning',false);
verify_core(root,here,{'monotone'}); abc_verify_adc;
if o.verify_only, fprintf('ABC_VERIFICATION_PASSED\n'); return; end
experiments=struct;
cachefile=outfile;
if o.resume && ~isfile(cachefile) && ~isempty(o.output_dir)
    cachefile=fullfile(here,'results','experiments_ABC.mat');
end
if o.resume && isfile(cachefile)
    saved=load(cachefile,'experiments','protocol');
    cached_protocol=saved.protocol;
    if isfield(cached_protocol,'B_increase') && ~isfield(cached_protocol,'B_increase_flags')
        assert(isequal(cached_protocol.B_increase,true),'Unrecognized prior B protocol.');
        cached_protocol=rmfield(cached_protocol,'B_increase');
        cached_protocol.B_increase_flags=[false true];
    end
    assert(abc_same_value(cached_protocol,protocol),'ABC protocol changed.');
    for key=o.phases
        if isfield(saved.experiments,key{1})
            experiments.(key{1})=saved.experiments.(key{1});
        end
    end
end
organization=struct('date','2026-10-16','source_cache',cachefile, ...
    'reused_actual_results',o.resume && isfile(cachefile), ...
    'note','Cached phase execution sources and timings are preserved.');
source_files={'main_line_search.m','TSITTCGA_LS.m'};
source_snapshot=cellfun(@(q)fileread(fullfile(here,q)),source_files,'UniformOutput',false);
matlab_version=version; timer=tic;
for phase=o.phases
    phase=phase{1}; assert(any(strcmp(phase,{'A','B','C'})));
    [cfg,solvers,params,labels,files]=abc_configs(phase,p,native,src);
    cache_valid=false;
    if isfield(experiments,phase) && ~strcmp(phase,'C')
        s=experiments.(phase);
        try
            abc_validate_cache(s,phase,p,meta,cases,cfg,params,labels,files,source_snapshot{2});
            cache_valid=true;
        catch cache_error
            warning('TSITTCGA:StaleCache', ...
                'Group %s will be recomputed because its cached inputs/source differ: %s', ...
                phase,cache_error.message);
        end
    end
    if cache_valid
        fprintf('ABC_REUSE %s\n',phase);
        if ~isempty(o.output_dir)
            elapsed=toc(timer);
            save(outfile,'experiments','protocol','organization','source_files','source_snapshot', ...
                'matlab_version','elapsed','o','-v7.3');
            abc_tables(experiments,outdir,'');
        end
        if o.plots, abc_plot(s,figdir); end
        continue
    end
    if strcmp(phase,'C')
        s=abc_import_c(src,cases,p,meta,cfg,params,labels,o.threads);
    else
        s=abc_run(phase,cases,p,meta,cfg,solvers,params,labels);
    end
    s.execution_sources=source_snapshot; s.ls_source=source_snapshot{2};
    s.solver_files=files; s.solver_sources=cellfun(@fileread,files,'UniformOutput',false);
    experiments.(phase)=s;
    elapsed=toc(timer);
    save(outfile,'experiments','protocol','source_files','source_snapshot', ...
        'matlab_version','elapsed','o','-v7.3'); % Preserve solved C before parity audit.
    if isfield(experiments,'B') && isfield(experiments,'C')
        b=experiments.B; c=experiments.C;
        assert(isequaln(b.NF(:,2),c.NF(:,1)) && isequaln(b.Itr(:,2),c.Itr(:,1)) ...
            && isequaln(b.RawNF(:,2),c.RawNF(:,1)) && isequaln(b.Success(:,2),c.Success(:,1)) ...
            && isequaln(b.Status(:,2),c.Status(:,1)) && isequaln(b.SolutionHash(:,2),c.SolutionHash(:,1)), ...
            'Full controlled core does not match native TSITTCGA.');
        experiments.B.native_full_core_parity=true; experiments.C.native_full_core_parity=true;
    end
    elapsed=toc(timer);
    save(outfile,'experiments','protocol','source_files','source_snapshot', ...
        'matlab_version','elapsed','o','-v7.3');
    if isempty(o.output_dir), abc_tables(experiments,outdir);
    else, abc_tables(experiments,outdir,''); end
    if o.plots, abc_plot(s,figdir); end
    fprintf('ABC_SAVED %s %.1fs\n',phase,toc(timer));
end
fprintf('ABC_STUDY_DONE %.1fs\n',toc(timer));
end

function abc_validate_cache(s,phase,p,meta,cases,cfg,params,labels,files,controlled_source)
% Compare the numerical experiment, not the layout of its orchestration code.
assert(strcmp(s.phase,phase) && abc_same_value(s.labels,labels), ...
    'Cached ABC methods or their order changed.');
assert(abc_same_value(s.params,p),'Cached TSITTCGA parameters changed.');
assert(strcmp(s.input_meta.original_solver,meta.original_solver), ...
    'The native TSITTCGA solver changed.');
assert(strcmp(s.input_meta.problem_source,meta.problem_source), ...
    'The problem mapping or its projection changed.');
assert(strcmp(s.ls_source,controlled_source),'Controlled core changed.');
assert(abc_same_value(s.cases,rmfield(cases,{'PR','x0'})), ...
    'The ABC case identities, dimensions, starts, or order changed.');
for field={'X0','problem_ids','dimensions'}
    key=field{1};
    assert(isfield(s.input_meta,key) && isfield(meta,key) ...
        && abc_same_value(s.input_meta.(key),meta.(key)), ...
        'The frozen ABC input field %s changed.',key);
end
assert(abc_same_value(s.configurations,cfg), ...
    'The ABC acceptance criteria or step-selection settings changed.');
assert(abc_same_value(s.method_params,params), ...
    'A complete C-study algorithm parameter changed.');
if strcmp(phase,'C')
    assert(numel(s.solver_sources)==numel(files), ...
        'The number of cached C solver sources changed.');
    for j=1:numel(files)
        current_source=fileread(files{j});
        cached_source=s.solver_sources{j};
        if strcmp(labels{j},'ADCGPM')
            current_source=abc_adc_source(current_source);
            cached_source=abc_adc_source(cached_source);
        end
        assert(strcmp(cached_source,current_source), ...
            'The C solver %s changed.',labels{j});
    end
end
% Execution snapshots and the original input provenance remain untouched.
end

function tf=abc_same_value(a,b)
% Reloaded anonymous handles retain their definition, not handle identity.
if isa(a,'function_handle') || isa(b,'function_handle')
    tf=isa(a,'function_handle') && isa(b,'function_handle') ...
        && strcmp(func2str(a),func2str(b));
elseif isstruct(a) && isstruct(b)
    fields=sort(fieldnames(a));
    tf=isequal(size(a),size(b)) && isequal(fields,sort(fieldnames(b)));
    if ~tf, return; end
    for i=1:numel(a)
        for j=1:numel(fields)
            key=fields{j};
            if ~abc_same_value(a(i).(key),b(i).(key)), tf=false; return; end
        end
    end
elseif iscell(a) && iscell(b)
    tf=isequal(size(a),size(b));
    if ~tf, return; end
    for i=1:numel(a)
        if ~abc_same_value(a{i},b{i}), tf=false; return; end
    end
else
    tf=isequaln(a,b);
end
end

function p=abc_main_parameters(main_source,iteration_limit)
% Only read the frozen parameter assignment; do not execute the driver.
block=regexp(main_source,'(?m)^\s*wjj\s*=\s*struct[\s\S]*?;','match','once');
assert(~isempty(block),'Cannot locate TSITTCGA parameters in the current main.m.');
o=struct('Itr_max',iteration_limit); %#ok<NASGU>
eval(block);
p=wjj;
end

function source=abc_adc_source(source)
% A local algorithm is independent of the surrounding experiment driver.
start=regexp(source,'(?m)^function[^\r\n]*=\s*ADCGPM\s*\(','start');
assert(isscalar(start),'Cannot locate exactly one ADCGPM function.');
source=source(start:end);
functions=regexp(source,'(?m)^function\s','start');
if numel(functions)>1, source=source(1:functions(2)-1); end
source=strrep(source,sprintf('\r\n'),sprintf('\n'));
source=regexprep(source,'\s+$','');
end

function [cfg,solvers,params,labels,files]=abc_configs(phase,p,native,src)
cfg={}; solvers={}; params={}; files={};
if any(strcmp(phase,{'A','B'}))
    base=struct('criterion','Proposed','sigma',p.sigma,'initial',p.rho_ls,'shrink',p.kappa, ...
        'warm',false,'increase',false,'lambda',.001,'nu',.8,'liu_mu',.01,'j_a',2,'j_l',.5,'label','Proposed');
    if strcmp(phase,'A')
        keys={'Proposed','Yin2021','Liu','Jiang'};
        labels={'Proposed (this paper)','Yin et al. [10]','Liu [15]','Jiang et al. [16]'};
        for j=1:4, q=base; q.criterion=keys{j}; q.label=labels{j}; cfg{j}=q; end
    else
        labels={'Fixed initial step','Warm start + increase'};
        for j=1:2, q=base; q.warm=j==2; q.increase=j==2; q.label=labels{j}; cfg{j}=q; end
    end
else
    labels={'TSITTCGA','AT2DFP','ADCGPM','TSIM-DFPM-PRP','HTTCGP'};
    solvers={@TSITTCGA,@AT2DFP,@ADCGPM,@TSIM_DFPM_PRP,@HTTCGP};
    adc=struct('Itr_max',p.Itr_max,'ell',.5,'sigma',.01,'shrink',.55, ...
        'rho1',10,'rho2',50,'mu',.4,'zeta',50,'xi',.01,'a',2,'gamma',1.9);
    for j=1:numel(labels)
        if j==1, params{j}=p;
        elseif j==3
            ix=find(strcmp(native.names,labels{j}));
            if isempty(ix)
                params{j}=adc;
            else
                assert(isscalar(ix) && abc_same_value(native.params{ix},adc), ...
                    'Formal ADCGPM parameters differ from the C-study configuration.');
                params{j}=native.params{ix};
            end
        else
            ix=find(strcmp(native.names,labels{j})); assert(isscalar(ix)); params{j}=native.params{ix};
        end
        if j==3
            file=fullfile(src,'ADCGPM.m');     % ADCGPM is a standalone solver file.
        else
            file=fullfile(src,[func2str(solvers{j}),'.m']);
        end
        assert(strcmp(functions(solvers{j}).file,file),'Wrong C solver path.');
        files{j}=file;
    end
end
end

function s=abc_import_c(src,cases,p,meta,cfg,params,labels,threads)
% Canonical C data are produced only by 单调方程组/main.m.
file=fullfile(src,'results','compare.mat'); d=load(file);
assert(d.matlab_threads==threads,'Rerun Group C with the requested thread count.');
assert(isequal(d.names,labels) && abc_same_value(d.params,params),'C parameters differ.');
assert(isfield(d,'SolutionHash'),'Rerun the native Group C main.m after the source update.');
for j=1:numel(d.source_files)
    assert(strcmp(fileread(fullfile(src,d.source_files{j})),d.source_snapshot{j}), ...
        'Canonical C results are stale. Rerun 单调方程组/main.m.');
end
nc=numel(cases); nm=numel(labels); M=nan(nc,nm);
s=struct('phase','C','params',p,'configurations',{cfg},'method_params',{params}, ...
    'labels',{labels},'cases',rmfield(cases,{'PR','x0'}),'input_meta',meta, ...
    'Itr',M,'NF',M,'Time',M,'LSNF',M,'NormF',M,'Feasibility',M, ...
    'Backtracks',M,'IncreaseTests',M,'IncreaseAccepts',M,'InitialAccepts',M, ...
    'Success',false(nc,nm),'Status',{cell(nc,nm)},'SolutionHash',{cell(nc,nm)}, ...
    'RawTime',M,'RawNF',M,'canonical_result',file,'wall',d.wall);
for i=1:nc
    q=find(strcmp({d.inst.name},cases(i).Problem) & [d.inst.n]==cases(i).n);
    r=cases(i).Start; assert(isscalar(q) && isequal(d.X0{q,r},cases(i).x0));
    for f={'Itr','NF','RawTime','RawNF','Success','Status','SolutionHash','Feasibility'}
        key=f{1}; s.(key)(i,:)=reshape(d.(key)(q,r,:),1,nm);
    end
    s.Time(i,:)=reshape(d.CPU(q,r,:),1,nm);
    s.NormF(i,:)=reshape(d.G(q,r,:),1,nm);
end
fprintf('ABC C imported from canonical native results.\n');
end

function s=abc_run(phase,cases,p,meta,cfg,solvers,params,labels)
nc=numel(cases); nm=numel(labels); M=nan(nc,nm);
s=struct('phase',phase,'params',p,'configurations',{cfg},'method_params',{params}, ...
    'labels',{labels},'cases',rmfield(cases,{'PR','x0'}),'input_meta',meta, ...
    'Itr',M,'NF',M,'Time',M,'LSNF',M,'NormF',M,'Feasibility',M, ...
    'Backtracks',M,'IncreaseTests',M,'IncreaseAccepts',M,'InitialAccepts',M, ...
    'Success',false(nc,nm),'Status',{cell(nc,nm)},'SolutionHash',{cell(nc,nm)}, ...
    'RawTime',M,'RawNF',M);
for j=1:nm
    short=p; if strcmp(phase,'C'), short=params{j}; end
    short.Itr_max=3;
    if strcmp(phase,'C'), solvers{j}(cases(1).PR,cases(1).x0,short);
    else, TSITTCGA_LS(cases(1).PR,cases(1).x0,short,cfg{j}); end
end
wall=tic;
for i=1:nc
    for j=mod(i-1+(0:nm-1),nm)+1
        if strcmp(phase,'C')
            [tc,nf,it,fn,d]=solvers{j}(cases(i).PR,cases(i).x0,params{j});
        else
            [tc,nf,it,fn,d]=TSITTCGA_LS(cases(i).PR,cases(i).x0,p,cfg{j});
            assert(d.actual_NF==1+d.inertia_nf+d.ls_nf+d.projection_nf);
            s.LSNF(i,j)=d.ls_nf; s.Backtracks(i,j)=d.backtracks;
            s.IncreaseTests(i,j)=d.increase_tests; s.IncreaseAccepts(i,j)=d.increase_accepts;
            s.InitialAccepts(i,j)=d.initial_accepts;
        end
        fx=cases(i).PR.f(d.x,1); normf=norm(fx);
        feas=norm(d.x-cases(i).PR.f(d.x,2));
        good=d.success && isreal(d.x) && all(isfinite(d.x)) && isreal(fx) ...
            && isfinite(normf) && normf<=1e-6 && isfinite(feas) && feas<=1e-6 ...
            && all(isfinite([tc,nf,it,fn])) && it<=p.Itr_max;
        assert(~d.success || good,'Invalid reported success.');
        s.Success(i,j)=good; s.Status{i,j}=d.status; s.NormF(i,j)=normf; s.Feasibility(i,j)=feas;
        s.RawTime(i,j)=d.actual_time; s.RawNF(i,j)=d.actual_NF;
        s.SolutionHash{i,j}=abc_hash(d.x);
        if good, s.NF(i,j)=nf; s.Itr(i,j)=it; s.Time(i,j)=tc; end
    end
    if mod(i,5)==0, fprintf('ABC %s %3d/%d %.1fs\n',phase,i,nc,toc(wall)); end
end
s.wall=toc(wall);
fprintf('ABC %s completed, unconditional means (NaN if incomplete)\n',phase);
for j=1:nm, fprintf('%-23s Itr=%.3f NF=%.3f Time=%.6f\n',labels{j},mean(s.Itr(:,j)),mean(s.NF(:,j)),mean(s.Time(:,j))); end
end

function h=abc_hash(x)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(double(x(:)),'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function abc_tables(experiments,outdir,prefix)
if nargin<3, prefix='abc_'; end
rows=cell(0,16); summary=cell(0,10); comparisons=cell(0,9); problemrows=cell(0,9);
for phase=fieldnames(experiments)'
    s=experiments.(phase{1}); nc=size(s.NF,1); nm=numel(s.labels);
    for j=1:nm
        for i=1:nc
            rows(end+1,:)={s.phase,s.cases(i).Problem,s.cases(i).n,s.cases(i).Start,s.labels{j}, ...
                s.Itr(i,j),s.NF(i,j),s.LSNF(i,j),s.Time(i,j),s.NormF(i,j),s.Feasibility(i,j), ...
                s.Backtracks(i,j),s.IncreaseTests(i,j),s.IncreaseAccepts(i,j),s.Success(i,j),s.Status{i,j}}; %#ok<AGROW>
        end
        ok=s.Success(:,j);
        summary(end+1,:)={s.phase,s.labels{j},nc,mean(ok),mean(s.Itr(:,j)),mean(s.NF(:,j)), ...
            mean(s.Time(:,j)),complete_mean(s.LSNF(:,j),ok),max(s.NormF(ok,j)),mean(s.NF(ok,j))}; %#ok<AGROW>
        ids=unique([[s.cases.n]',str2double(erase(string({s.cases.Problem}'),'P'))],'rows','stable');
        for k=1:size(ids,1)
            mask=[s.cases.n]'==ids(k,1) & strcmp({s.cases.Problem}',sprintf('P%d',ids(k,2)));
            problemrows(end+1,:)={s.phase,sprintf('P%d',ids(k,2)),ids(k,1),s.labels{j},sum(mask), ...
                mean(s.Itr(mask,j)),mean(s.NF(mask,j)),mean(s.Time(mask,j)),max(s.NormF(mask,j))}; %#ok<AGROW>
        end
    end
    pairs=[ones(nm-1,1),(2:nm)']; if strcmp(s.phase,'B'), pairs=[2 1]; end
    for k=1:size(pairs,1)
        a=pairs(k,1); b=pairs(k,2); mask=s.Success(:,a)&s.Success(:,b);
        comparisons(end+1,:)={s.phase,[s.labels{a},' / ',s.labels{b}],sum(mask), ...
            paired_ratio(s.NF(mask,a),s.NF(mask,b)),paired_ratio(s.Itr(mask,a),s.Itr(mask,b)), ...
            paired_ratio(s.Time(mask,a),s.Time(mask,b)),paired_ratio(s.LSNF(mask,a),s.LSNF(mask,b)), ...
            mean(s.Success(:,a)),mean(s.Success(:,b))}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'Group','Problem','n','Start','Method','Itr','NF','LSNF', ...
    'Time','NormF','Feasibility','Backtracks','IncreaseTests','IncreaseAccepts','Success','Status'});
writetable(T,fullfile(outdir,[prefix,'all_runs.csv']),'Encoding','UTF-8');
T=cell2table(summary,'VariableNames',{'Group','Method','Instances','SuccessRate','MeanItr','MeanNF', ...
    'MeanTime','MeanLSNF','MaxNormFSolved','MeanNFSolved'});
writetable(T,fullfile(outdir,[prefix,'summary.csv']),'Encoding','UTF-8');
T=cell2table(comparisons,'VariableNames',{'Group','Comparison','CommonInstances','NF_Ratio','Itr_Ratio', ...
    'Time_Ratio','LSNF_Ratio','SuccessRateA','SuccessRateB'});
writetable(T,fullfile(outdir,[prefix,'paired_comparisons.csv']),'Encoding','UTF-8');
T=cell2table(problemrows,'VariableNames',{'Group','Problem','n','Method','Starts','MeanItr','MeanNF','MeanTime','MaxNormF'});
writetable(T,fullfile(outdir,[prefix,'problem_summary.csv']),'Encoding','UTF-8');
end

function abc_plots(experiments,folder)
for key=fieldnames(experiments)', abc_plot(experiments.(key{1}),folder); end
end

function abc_plot(s,folder)
% One profile case is one problem--dimension pair, averaged over five starts.
% A pair with any failed attempt receives infinite cost for that method.
% The panels are drawn by profile_panel.m, the same routine used by the
% native Group C figure, so every profile figure of the paper shares one
% style block (line width, markers, truncation at rho=1, no overlap).
labels=s.labels;
labels=regexprep(labels,'^Proposed$','Proposed (this paper)');
labels=regexprep(labels,'^Yin \(2021\)$','Yin [2021]');
labels=regexprep(labels,'^Liu \(2026\)$','Liu [2026]');
labels=regexprep(labels,'^Jiang et al\. \(2026\)$','Jiang [2026]');
% Manuscript numbering of the cited sources (order of the reference list).
labels=regexprep(labels,'^Yin \[2021\]$','Yin et al. [10]');
labels=regexprep(labels,'^Liu \[2026\]$','Liu [15]');
labels=regexprep(labels,'^Jiang \[2026\]$','Jiang et al. [16]');
addpath(fileparts(fileparts(mfilename('fullpath'))),'-end');   % profile_panel.m
problem=string({s.cases.Problem}'); dimensions=[s.cases.n]';
keys=unique(table(problem,dimensions),'rows','stable');
np=height(keys); nm=numel(labels); fields={'Itr','NF','Time'};
means=struct('Itr',nan(np,nm),'NF',nan(np,nm),'Time',nan(np,nm));
success=false(np,nm);
for i=1:np
    mask=problem==keys.problem(i) & dimensions==keys.dimensions(i);
    starts=sort([s.cases(mask).Start]);
    assert(isequal(starts,1:5),'Expected exactly five distinct starts per problem-dimension pair.');
    success(i,:)=all(s.Success(mask,:),1);
    for m=1:numel(fields)
        runs=s.(fields{m})(mask,:); ok=s.Success(mask,:);
        assert(all(isfinite(runs(ok))) && all(runs(ok)>0), ...
            'A successful run has an invalid %s cost.',fields{m});
        means.(fields{m})(i,:)=mean(runs,1);
    end
end
fprintf('Group %s profiles: %d problem-dimension pairs, mean over 5 starts per pair.\n',s.phase,np);
fprintf('Group %s complete-pair success fractions %s\n',s.phase,mat2str(mean(success,1),4));
f=profile_figure();
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
fields={'Itr','NF','Time'}; titles={'(a) Itr','(b) NF','(c) Tcpu'};
for m=1:3
    ax=nexttile; hold(ax,'on');
    cost=means.(fields{m}); cost(~success | ~isfinite(cost) | cost<=0)=Inf;
    best=min(cost,[],2); ratios=cost./best; ratios(~isfinite(ratios))=Inf;
    h=profile_panel(ax,ratios,labels,titles{m},m==1);
    legend(ax,'off');
    if m==1, firstHandle=h; end
    fprintf('Group %s %s best fractions %s\n',s.phase,fields{m},mat2str(mean(abs(ratios-1)<1e-12,1),4));
end
figure_frame(tl,firstHandle,labels,sprintf('Group %s: five projection algorithms, %d problem-dimension pairs',s.phase,np),5);
exportOptions={'Units','centimeters','Width',f.Position(3),'Height',f.Position(4), ...
    'Padding','figure','PreserveAspectRatio','on'};
exportgraphics(f,fullfile(folder,['group_',s.phase,'_profiles.pdf']),'ContentType','vector',exportOptions{:});
exportgraphics(f,fullfile(folder,['group_',s.phase,'_profiles.png']),'Resolution',300,exportOptions{:});
end

function f=profile_figure()
% Canvas for every 1x3 Dolan-More profile (pbaspect 1.2 on each panel).
% The panels are limited by the tile WIDTH, not by the canvas height, so a
% canvas taller than about 5.8 cm only adds white bands above and below the
% plot; 18 x 5.8 cm leaves ~10 bp of margin instead of the ~98 bp that
% 18 x 7.5 cm produced. Change the size here only -- export_figure style
% calls read Width/Height from f.Position, so the two never drift apart.
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 5.3]);
end

function abc_restore(savedpath,savedrng,savedthreads)
path(savedpath); rng(savedrng); maxNumCompThreads(savedthreads);
end
function abc_verify_adc
% Independent first-step calculation for F(x)=2x-1 on the nonnegative orthant.
pr=struct('name','AffineADCGPMCheck','f',@abc_affine_problem);
q=struct('Itr_max',1);
[~,~,~,~,d]=ADCGPM(pr,ones(3,1),q);
assert(norm(d.x-.42525*ones(3,1))<1e-14 && d.actual_NF==5 ...
    && d.accepted_updates==1 && d.max_backtracks==2,'ADCGPM first-step equation mismatch.');
q.Itr_max=2000;
[~,nf,it,fn,d]=ADCGPM(pr,ones(3,1),q);
assert(d.success && fn<=1e-6 && norm(d.x-.5*ones(3,1))<=1e-6 ...
    && nf==d.actual_NF && it<=2000,'ADCGPM affine convergence check failed.');
fprintf('ADCGPM first-step projection and affine convergence verified.\n');
end

function out=abc_affine_problem(x,mode)
if mode==1, out=2*x-1; else, out=max(x,0); end
end
function fair_study(options)
% Independent A-study: common-sigma sensitivity, equal tuning budgets,
% then one frozen configuration per criterion on 200 fresh test starts.
o=struct('resume',true,'figs_only',false,'verify_only',false,'plots',true,'threads',1);
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'Unknown fair-study option %s.',key{1}); o.(key{1})=options.(key{1});
end
validateattributes(o.threads,{'numeric'},{'scalar','integer','positive','finite'});
here=fileparts(mfilename('fullpath')); root=fileparts(here);
outdir=fullfile(here,'results'); figdir=fullfile(here,'figures');
outfile=fullfile(outdir,'A_fair.mat');
if o.figs_only
    d=load(outfile,'fair');
    assert(d.fair.complete,'Finish the frozen study before plotting final profiles.');
    fair_plots(d.fair,figdir); return
end
oldpath=path; oldrng=rng; oldthreads=maxNumCompThreads;
cleanup=onCleanup(@()abc_restore(oldpath,oldrng,oldthreads)); %#ok<NASGU>
maxNumCompThreads(o.threads);
src=source_dir(root,'monotone'); addpath(src,'-begin'); addpath(here,'-begin');
assert(strcmp(functions(@TSITTCGA_LS).file,fullfile(here,'TSITTCGA_LS.m')), ...
    'The fair study must use the controlled solver in this directory.');
clear problems; pf=str2func('problems');
assert(strcmp(functions(pf).file,fullfile(src,'problems.m')),'Wrong monotone problems source.');
historical=load(fullfile(src,'results','compare.mat'),'wjj','inst','X0'); p=historical.wjj;
assert(abc_same_value(p,abc_main_parameters(fileread(fullfile(src,'main.m')),p.Itr_max)), ...
    'Frozen TSITTCGA parameters do not match the current monotone main.m.');
assert(p.Itr_max==2000 && numel(historical.inst)==40 && size(historical.X0,2)==5, ...
    'The fair study requires the original eight-problem, five-dimension design.');
protocol=fair_protocol(o.threads);
cfg=fair_candidates(p,protocol); fair_verify(p,cfg);
if o.verify_only, fprintf('FAIR_A_VERIFICATION_PASSED\n'); return; end
sources=struct('controlled',fileread(fullfile(here,'TSITTCGA_LS.m')), ...
    'problems',fileread(fullfile(src,'problems.m')), ...
    'native',fileread(fullfile(src,'TSITTCGA.m')));
history_hashes=cellfun(@abc_hash,historical.X0,'UniformOutput',false);
if isfile(outfile)
    assert(o.resume || o.figs_only, ...
        'A_fair.mat already exists. Use resume=true to preserve the frozen study.');
    d=load(outfile,'fair'); fair=d.fair;
    assert(abc_same_value(fair.protocol,protocol),'The frozen fair-study protocol changed.');
    assert(abc_same_value(fair.base_params,p),'The frozen non-line-search parameters changed.');
    assert(abc_same_value(fair.candidates,cfg),'The frozen line-search candidate set changed.');
    assert(abc_same_value(fair.sources,sources),'A controlled solver or problem source changed.');
    assert(abc_same_value(fair.original_instances,historical.inst) ...
        && isequal(fair.historical_start_hashes,history_hashes), ...
        'The historical reference inputs changed.');
    assert(strcmp(fair.input_digest,fair_input_digest(fair.inputs)), ...
        'A frozen development or test input changed.');
    if fair.selection_frozen
        [selected,scores]=fair_select(fair.development,protocol);
        assert(isequal(selected,fair.selected_idx) ...
            && isequaln(scores,fair.selection_scores), ...
            'Frozen selected configurations no longer match the development data.');
        for j=1:4
            expected=cfg{selected(j),j}; [effective,~]=fair_effective(p,expected);
            assert(abc_same_value(fair.selected_cfg{j},expected) ...
                && abc_same_value(fair.effective_params{j},effective), ...
                'The actually selected line-search configuration differs from its frozen index.');
        end
    end
else
    assert(~o.figs_only,'There are no fair-study data to plot.');
    fair=struct('protocol',protocol,'base_params',p,'candidates',{cfg}, ...
        'sources',sources,'execution_driver',fileread([mfilename('fullpath'),'.m']), ...
        'matlab_version',version,'historical_start_hashes',{history_hashes}, ...
        'original_instances',historical.inst,'selection_frozen',false, ...
        'selected_idx',[],'selected_cfg',{{}},'effective_params',{{}}, ...
        'selection_scores',[],'complete',false);
    fair.inputs=fair_inputs(historical,pf,protocol,history_hashes);
    fair.input_digest=fair_input_digest(fair.inputs);
    fair.development=fair_empty(32,24,4,1);
    fair.test=fair_empty(200,1,4,3);
    fair.created_at=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss z'));
    if ~isfolder(outdir), mkdir(outdir); end
    save(outfile,'fair','-v7.3','-nocompression'); % Freeze inputs and candidates before any performance observation.
    fair_tables(fair,outdir);
    fprintf('FAIR_A_FROZEN 32 development / 200 fresh test cases; 24 candidates per criterion.\n');
end
clear historical;
if o.figs_only
    fair_tables(fair,outdir); if o.plots, fair_plots(fair,figdir); end; return
end
if fair.complete
    fprintf('FAIR_A_REUSE_COMPLETE\n'); fair_tables(fair,outdir);
    if o.plots, fair_plots(fair,figdir); end; return
end
% Equal untimed warm-up for all four criteria; no pilot result selects parameters.
warm=fair_problem(pf,1,9); short=p; short.Itr_max=3;
for j=1:4, TSITTCGA_LS(warm,ones(9,1),short,cfg{1,j}); end
wall=tic;
for candidate=1:24
    for j=mod(candidate-1+(0:3),4)+1
        if fair.development.ConfigDone(candidate,j), continue; end
        [effective,ls]=fair_effective(p,cfg{candidate,j});
        for i=1:32
            item=fair.inputs.development.cases(i);
            pr=fair_problem(pf,item.id,item.n);
            x0=fair.inputs.development.X0{item.InputRow,item.Start};
            record=fair_run(pr,x0,effective,ls);
            fair.development=fair_record(fair.development,i,candidate,j,1,record);
        end
        fair.development.ConfigDone(candidate,j)=true;
        fair.updated_at=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss z'));
        save(outfile,'fair','-v7.3','-nocompression');
        fprintf('FAIR_A_DEV candidate %2d/24 criterion %d/4, success %d/32, elapsed %.1fs\n', ...
            candidate,j,sum(fair.development.Success(:,candidate,j)),toc(wall));
    end
    if candidate==5
        fair_tables(fair,outdir);
        if o.plots, fair_plots(fair,figdir); end
        fprintf('FAIR_A_COMMON_SIGMA_DONE (the same development runs remain in the 24-candidate budget).\n');
    end
end
assert(all(fair.development.ConfigDone,'all'),'Development budget is incomplete.');
if ~fair.selection_frozen
    [fair.selected_idx,fair.selection_scores]=fair_select(fair.development,protocol);
    for j=1:4
        fair.selected_cfg{j}=cfg{fair.selected_idx(j),j};
        [fair.effective_params{j},~]=fair_effective(p,fair.selected_cfg{j});
    end
    fair.selection_frozen=true;
    fair.selection_frozen_at=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss z'));
    save(outfile,'fair','-v7.3','-nocompression'); % No formal-test run is permitted before this saved freeze.
    fair_tables(fair,outdir);
    fprintf('FAIR_A_SELECTION_FROZEN indices [%s]\n',num2str(fair.selected_idx));
end
for j=1:4
    q=fair.effective_params{j}; q.Itr_max=3;
    TSITTCGA_LS(warm,ones(9,1),q,fair.selected_cfg{j});
end
for i=1:200
    item=fair.inputs.test.cases(i); pr=fair_problem(pf,item.id,item.n);
    x0=fair.inputs.test.X0{item.InputRow,item.Start};
    for rep=1:3
        for j=mod(i+rep-2+(0:3),4)+1
            if fair.test.Done(i,1,j), continue; end
            rec=fair_run(pr,x0,fair.effective_params{j},fair.selected_cfg{j});
            if rep>1, fair_repeat_check(fair.test,i,j,rec); end
            fair.test=fair_record(fair.test,i,1,j,rep,rec);
            if rep==3
                fair.test.Done(i,1,j)=true;
                if fair.test.Success(i,1,j)
                    fair.test.Time(i,1,j)=median(fair.test.RawTime(i,1,j,:));
                end
            end
        end
    end
    if mod(i,5)==0
        fair.updated_at=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss z'));
        save(outfile,'fair','-v7.3','-nocompression');
        fprintf('FAIR_A_TEST %3d/200 cases (three repetitions), elapsed %.1fs\n',i,toc(wall));
    end
end
assert(all(fair.test.Done,'all'),'Formal test is incomplete.');
fair.complete=true; fair.test.ConfigDone(:)=true;
fair.completed_at=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss z'));
save(outfile,'fair','-v7.3','-nocompression'); fair_tables(fair,outdir);
if o.plots, fair_plots(fair,figdir); end
fprintf('FAIR_A_STUDY_DONE %.1fs\n',toc(wall));
end

function protocol=fair_protocol(threads)
protocol=struct('version',1,'date','2026-10-15','threads',threads, ...
    'keys',{{'Proposed','Yin2021','Liu','Jiang'}}, ...
    'labels',{{'Proposed','Yin et al. (2020)','Liu (2026)','Jiang et al. (2026)'}}, ...
    'development_seed',83117,'test_seed',94421,'candidate_seed',72563, ...
    'development_dimensions',[1000 3000],'development_dimensions6',[200 500], ...
    'test_dimensions',[10000 30000 50000 80000 100000], ...
    'test_dimensions6',[1000 3000 5000 8000 10000], ...
    'development_starts',2,'test_starts',5,'candidate_count',24, ...
    'sigma_grid',[1e-4 1e-2 .1 .5 .8],'initial',.8,'shrink',.5, ...
    'warm',false,'increase',false,'repetitions',3,'iteration_limit',2000, ...
    'epsilon',1e-6,'epsilon_C',1e-6,'selection_failure_cost',1e8, ...
    'sigma_range',[1e-4 .9],'proposed_mu_range',[.05 1], ...
    'proposed_a_range',[1e-4 1],'yin_lambda_range',[1e-4 .1], ...
    'yin_nu_range',[.1 2],'liu_mu_range',[1e-4 1], ...
    'jiang_ell_range',[.25 1.5],'jiang_a_range',[1.1 10], ...
    'selection_rule','maximum success count; minimum all-case mean log(NF or 1e8); smallest candidate index', ...
    'test_time','median of three actual solver times; rotating order', ...
    'development_role','select parameters and common-sigma sensitivity only', ...
    'test_role','fresh starts, held out until configuration freeze');
end

function cfg=fair_candidates(p,protocol)
saved=rng; cleanup=onCleanup(@()rng(saved)); %#ok<NASGU>
rng(protocol.candidate_seed,'twister');
random_sigma=exp(log(1e-4)+rand(19,1)*(log(.9)-log(1e-4)));
base=struct('criterion','Proposed','sigma',.8,'initial',.8,'shrink',.5, ...
    'warm',false,'increase',false,'lambda',.001,'nu',.8,'liu_mu',.01, ...
    'j_a',2,'j_l',.5,'label','Proposed','proposed_mu',p.mu,'proposed_a',p.a);
cfg=cell(24,4);
for j=1:4
    for c=1:24
        ls=base; ls.criterion=protocol.keys{j}; ls.label=protocol.labels{j};
        if c<=5
            ls.sigma=protocol.sigma_grid(c);
        else
            ls.sigma=random_sigma(c-5);
            switch j
                case 1
                    ls.proposed_mu=.05+.95*rand;
                    ls.proposed_a=exp(log(1e-4)+rand*log(1e4));
                case 2
                    ls.lambda=exp(log(1e-4)+rand*log(1e3));
                    ls.nu=exp(log(.1)+rand*log(20));
                case 3, ls.liu_mu=exp(log(1e-4)+rand*log(1e4));
                case 4
                    ls.j_l=.25+1.25*rand;
                    ls.j_a=exp(log(1.1)+rand*(log(10)-log(1.1)));
            end
        end
        cfg{c,j}=ls;
        fair_effective(p,ls);
    end
end
end

function [effective,ls]=fair_effective(p,configuration)
effective=p; ls=configuration;
if strcmp(ls.criterion,'Proposed')
    effective.mu=ls.proposed_mu; effective.a=ls.proposed_a;
end
effective.sigma=ls.sigma;
assert(abc_same_value(rmfield(effective,{'sigma','mu','a'}),rmfield(p,{'sigma','mu','a'})), ...
    'A non-line-search parameter was altered.');
if ~strcmp(ls.criterion,'Proposed')
    assert(effective.mu==p.mu && effective.a==p.a,'Unused proposed-criterion parameters changed.');
end
assert(ls.sigma>0 && ls.sigma<1 && ls.initial==.8 && ls.shrink==.5 ...
    && ~ls.warm && ~ls.increase && ls.lambda>0 && ls.lambda<=ls.nu ...
    && ls.liu_mu>0 && ls.j_a>1 && ls.j_l>0 ...
    && effective.mu>0 && effective.mu<=1 && effective.a>0, ...
    'Invalid or noncomparable line-search configuration.');
end

function fair_verify(p,cfg)
assert(isequal(size(cfg),[24 4]),'Every criterion requires exactly 24 candidates.');
for j=1:4
    textkeys=cell(24,1);
    for c=1:24
        [q,ls]=fair_effective(p,cfg{c,j});
        assert(ls.sigma>=1e-4 && ls.sigma<=.9 && ls.proposed_mu>=.05 && ls.proposed_mu<=1 ...
            && ls.proposed_a>=1e-4 && ls.proposed_a<=1 ...
            && ls.lambda>=1e-4 && ls.lambda<=.1 && ls.nu>=.1 && ls.nu<=2 ...
            && ls.liu_mu>=1e-4 && ls.liu_mu<=1 && ls.j_l>=.25 && ls.j_l<=1.5 ...
            && ls.j_a>=1.1 && ls.j_a<=10,'A candidate is outside the frozen protocol ranges.');
        textkeys{c}=jsonencode(ls);
        assert(q.rho==p.rho && q.varsigma==p.varsigma && q.eta==p.eta ...
            && q.tbar==p.tbar && q.c==p.c && q.r==p.r,'Direction/inertia changed.');
    end
    assert(numel(unique(textkeys))==24,'Duplicate candidate configurations.');
end
for c=1:24
    assert(numel(unique(cellfun(@(s)s.sigma,cfg(c,:))))==1,'The common sigma samples differ.');
end
pr=struct('name','FairAffineCheck','f',@abc_affine_problem);
for j=1:4
    [q,ls]=fair_effective(p,cfg{2,j});
    r=fair_run(pr,ones(3,1),q,ls);
    assert(r.Success && r.NormF<=1e-6 && r.Feasibility<=1e-6,'Affine verification failed.');
    [~,nf,~,~,d]=TSITTCGA_LS(pr,.5*ones(3,1),q,ls);
    assert(d.success && nf==1 && d.actual_NF==1,'Initial-root or NF safeguard failed.');
end
% Exercise new input generation, projection, checkpoint masks and selection
% using tiny independent fixtures. These checks do not write experiment files.
protocol=fair_protocol(1); fixture=protocol;
fixture.development_dimensions=[9 11]; fixture.development_dimensions6=[13 17];
fixture.test_dimensions=[19 23 29 31 37]; fixture.test_dimensions6=[41 43 47 53 59];
historical=struct('inst',struct('id',{},'n',{}));
for id=1:8
    dims=fixture.test_dimensions; if id==6, dims=fixture.test_dimensions6; end
    for n=dims, historical.inst(end+1)=struct('id',id,'n',n); end %#ok<AGROW>
end
inputs=fair_inputs(historical,@problems,fixture,{});
assert(ischar(fair_input_digest(inputs)) && numel(inputs.test.cases)==200);
s=fair_empty(32,24,4,1); s.ConfigDone(:)=true; s.Done(:)=true;
s.Success(:)=true; s.NF(:)=10;
s.Success(1,1,:)=false; s.NF(:,2,:)=5; s.NF(:,3,:)=5;
[selected,~]=fair_select(s,protocol);
assert(isequal(selected,[2 2 2 2]),'Selection must maximize success, minimize all-case cost, then use index.');
s=fair_empty(1,1,4,3);
[q,ls]=fair_effective(p,cfg{2,1}); rec=fair_run(pr,ones(3,1),q,ls);
s=fair_record(s,1,1,1,1,rec);
assert(~s.Done(1,1,1),'One timing repetition is not a completed formal run.');
fair_repeat_check(s,1,1,rec);
s=fair_record(s,1,1,1,2,rec); s=fair_record(s,1,1,1,3,rec);
assert(s.RepeatVerified(1,1,1),'Timing repetitions were not verified.');
fprintf('Fair-study configurations, affine checks, NF accounting, fresh inputs, masks and selection checked.\n');
end

function inputs=fair_inputs(historical,pf,protocol,history_hashes)
saved=rng; cleanup=onCleanup(@()rng(saved)); %#ok<NASGU>
inputs=struct;
all_new_hashes=cell(0,1);
for group={'development','test'}
    name=group{1}; isdev=strcmp(name,'development');
    if isdev
        rng(protocol.development_seed,'twister'); nr=2;
    else
        rng(protocol.test_seed,'twister'); nr=5;
    end
    entries=struct('id',{},'n',{});
    for id=1:8
        if isdev
            dims=protocol.development_dimensions;
            if id==6, dims=protocol.development_dimensions6; end
        else
            dims=protocol.test_dimensions;
            if id==6, dims=protocol.test_dimensions6; end
            ids=[historical.inst.id]; actual=[historical.inst(ids==id).n];
            assert(isequal(actual(:)',dims),'Unexpected original problem dimensions.');
        end
        for n=dims, entries(end+1)=struct('id',id,'n',n); end %#ok<AGROW>
    end
    X0=cell(numel(entries),nr); hashes=X0;
    cases=struct('id',{},'n',{},'Start',{},'InputRow',{});
    for row=1:numel(entries)
        pr=fair_problem(pf,entries(row).id,entries(row).n);
        for r=1:nr
            if entries(row).id==3, x=-1+2*rand(entries(row).n,1);
            else, x=4*rand(entries(row).n,1); end
            x=pr.f(x,2);
            assert(all(isfinite(x)) && norm(x-pr.f(x,2))<=1e-6,'Invalid generated feasible start.');
            h=abc_hash(x);
            assert(~any(strcmp(h,history_hashes(:))) && ~any(strcmp(h,all_new_hashes)), ...
                'Development/test starts must be new and distinct.');
            X0{row,r}=x; hashes{row,r}=h; all_new_hashes{end+1,1}=h; %#ok<AGROW>
            cases(end+1)=struct('id',entries(row).id,'n',entries(row).n,'Start',r,'InputRow',row); %#ok<AGROW>
        end
    end
    inputs.(name)=struct('entries',entries,'X0',{X0},'hashes',{hashes},'cases',cases);
end
assert(numel(inputs.development.cases)==32 && numel(inputs.test.cases)==200);
end

function digest=fair_input_digest(inputs)
hashes=cell(0,1);
for group={'development','test'}
    q=inputs.(group{1});
    for i=1:numel(q.X0)
        h=abc_hash(q.X0{i});
        assert(strcmp(h,q.hashes{i}),'A stored start does not match its original SHA256.');
        hashes{end+1,1}=h; %#ok<AGROW>
    end
    hashes{end+1,1}=jsonencode(q.cases); %#ok<AGROW>
end
digest=abc_hash(double(unicode2native(strjoin(hashes,'|'),'UTF-8')));
end

function pr=fair_problem(pf,id,n)
pr=struct('name',sprintf('P%d',id),'n',n,'f',@(x,mode)pf(x,mode,id));
end

function s=fair_empty(nc,ns,nm,reps)
M=nan(nc,ns,nm);
s=struct('Itr',M,'NF',M,'Time',M,'NormF',M,'Feasibility',M,'LSNF',M, ...
    'Success',false(nc,ns,nm),'Done',false(nc,ns,nm), ...
    'Status',{cell(nc,ns,nm)},'SolutionHash',{cell(nc,ns,nm)}, ...
    'RawNF',nan(nc,ns,nm,reps),'RawTime',nan(nc,ns,nm,reps), ...
    'RawItr',nan(nc,ns,nm,reps),'RepeatVerified',false(nc,ns,nm), ...
    'ConfigDone',false(ns,nm));
end

function r=fair_run(pr,x0,p,ls)
[tc,nf,it,fn,d]=TSITTCGA_LS(pr,x0,p,ls);
assert(d.actual_NF==1+d.inertia_nf+d.ls_nf+d.projection_nf,'Inconsistent actual NF.');
fx=pr.f(d.x,1); normf=norm(fx); feas=norm(d.x-pr.f(d.x,2));
good=d.success && isreal(d.x) && all(isfinite(d.x)) && isreal(fx) ...
    && isfinite(normf) && normf<=1e-6 && isfinite(feas) && feas<=1e-6 ...
    && all(isfinite([tc,nf,it,fn])) && it<=p.Itr_max;
assert(~d.success || good,'A reported success violates the independent stopping check.');
if good
    assert(nf==d.actual_NF && abs(fn-normf)<=1e-12*max(1,normf),'Invalid successful output or NF count.');
else
    nf=NaN; it=NaN; tc=NaN;
end
r=struct('Success',good,'NF',nf,'Itr',it,'Time',tc,'NormF',normf, ...
    'Feasibility',feas,'LSNF',d.ls_nf,'Status',d.status, ...
    'SolutionHash',abc_hash(d.x),'RawNF',d.actual_NF,'RawTime',d.actual_time, ...
    'RawItr',d.loop_iterations);
end

function s=fair_record(s,i,c,j,rep,r)
if rep==1
    for field={'Success','NF','Itr','Time','NormF','Feasibility','LSNF'}
        key=field{1}; s.(key)(i,c,j)=r.(key);
    end
    s.Status{i,c,j}=r.Status; s.SolutionHash{i,c,j}=r.SolutionHash;
end
s.RawNF(i,c,j,rep)=r.RawNF; s.RawTime(i,c,j,rep)=r.RawTime; s.RawItr(i,c,j,rep)=r.RawItr;
if size(s.RawNF,4)==1, s.Done(i,c,j)=true; end
if rep==3, s.RepeatVerified(i,c,j)=true; end
end

function fair_repeat_check(s,i,j,r)
assert(isequaln(s.Success(i,1,j),r.Success) && strcmp(s.Status{i,1,j},r.Status) ...
    && isequaln(s.NF(i,1,j),r.NF) && isequaln(s.Itr(i,1,j),r.Itr) ...
    && s.RawNF(i,1,j,1)==r.RawNF && s.RawItr(i,1,j,1)==r.RawItr ...
    && strcmp(s.SolutionHash{i,1,j},r.SolutionHash), ...
    'Deterministic NF, Itr, status, success or returned-point SHA256 changed across timing repetitions.');
end

function [selected,scores]=fair_select(s,protocol)
assert(all(s.ConfigDone,'all') && all(s.Done,'all'),'Cannot select from an incomplete budget.');
selected=zeros(1,4); scores=nan(24,4,2);
for j=1:4
    ranking=nan(24,3);
    for c=1:24
        success=s.Success(:,c,j); costs=repmat(protocol.selection_failure_cost,32,1);
        costs(success)=s.NF(success,c,j);
        count=sum(success); logcost=mean(log(costs));
        scores(c,j,1)=count; scores(c,j,2)=logcost;
        ranking(c,:)=[-count logcost c];
    end
    ranking=sortrows(ranking,[1 2 3]); selected(j)=ranking(1,3);
end
end

function fair_tables(fair,outdir)
labels=fair.protocol.labels; cfg=fair.candidates;
candidates=cell(0,15); development=cell(0,13); test=cell(0,13); summary=cell(0,12); pairs=cell(0,8);
for j=1:4
    for c=1:24
        q=cfg{c,j}; selected=fair.selection_frozen && fair.selected_idx(j)==c;
        completed=fair.development.ConfigDone(c,j); n=NaN; score=NaN;
        if completed
            ok=fair.development.Success(:,c,j); costs=repmat(fair.protocol.selection_failure_cost,32,1);
            costs(ok)=fair.development.NF(ok,c,j); n=sum(ok); score=mean(log(costs));
        end
        candidates(end+1,:)={labels{j},c,q.sigma,q.proposed_mu,q.proposed_a,q.lambda,q.nu, ...
            q.liu_mu,q.j_l,q.j_a,q.initial,q.shrink,completed,n,selected}; %#ok<AGROW>
        for i=1:32
            item=fair.inputs.development.cases(i); s=fair.development;
            development(end+1,:)={labels{j},c,item.id,item.n,item.Start,s.Done(i,c,j), ...
                s.Success(i,c,j),s.Itr(i,c,j),s.NF(i,c,j),s.RawNF(i,c,j,1), ...
                s.NormF(i,c,j),s.Status{i,c,j},score}; %#ok<AGROW>
        end
    end
    s=fair.test;
    for i=1:200
        item=fair.inputs.test.cases(i);
        test(end+1,:)={labels{j},item.id,item.n,item.Start,s.Done(i,1,j),s.Success(i,1,j), ...
            s.Itr(i,1,j),s.NF(i,1,j),s.Time(i,1,j),s.NormF(i,1,j),s.Feasibility(i,1,j), ...
            s.Status{i,1,j},s.SolutionHash{i,1,j}}; %#ok<AGROW>
    end
    done=s.Done(:,1,j); ok=s.Success(:,1,j); solved=done&ok;
    meanitr=NaN; meannf=NaN; meantime=NaN; successrate=NaN; selected=NaN; maxfn=NaN; maxfeas=NaN;
    if fair.selection_frozen, selected=fair.selected_idx(j); end
    if any(solved), maxfn=max(s.NormF(solved,1,j)); maxfeas=max(s.Feasibility(solved,1,j)); end
    if all(done)
        successrate=mean(ok);
        if all(ok), meanitr=mean(s.Itr(:,1,j)); meannf=mean(s.NF(:,1,j)); meantime=mean(s.Time(:,1,j)); end
    end
    summary(end+1,:)={labels{j},selected,sum(done),200,successrate,meanitr,meannf,meantime, ...
        maxfn,maxfeas,all(s.RepeatVerified(:,1,j)),fair.complete}; %#ok<AGROW>
    if j>1 && all(s.Done(:,1,1)) && all(done)
        common=s.Success(:,1,1)&ok;
        pairs(end+1,:)={labels{1},labels{j},sum(common),200, ...
            fair_georatio(s.NF(common,1,1),s.NF(common,1,j)), ...
            fair_georatio(s.Itr(common,1,1),s.Itr(common,1,j)), ...
            fair_georatio(s.Time(common,1,1),s.Time(common,1,j)), ...
            'ratios use only paired common successes; profiles use all 200'}; %#ok<AGROW>
    end
end
writetable(cell2table(candidates,'VariableNames',{'Method','Candidate','Sigma','ProposedMu','ProposedA', ...
    'YinLambda','YinNu','LiuMu','JiangEll','JiangA','Initial','Shrink','DevelopmentComplete','DevelopmentSuccesses','Selected'}), ...
    fullfile(outdir,'A_fair_candidates.csv'),'Encoding','UTF-8');
writetable(cell2table(development,'VariableNames',{'Method','Candidate','Problem','n','Start','Done','Success', ...
    'Itr','NF','RawNF','NormF','Status','AllCaseMeanLogCost'}),fullfile(outdir,'A_fair_development.csv'),'Encoding','UTF-8');
writetable(cell2table(test,'VariableNames',{'Method','Problem','n','Start','Done','Success','Itr','NF','MedianTime', ...
    'NormF','Feasibility','Status','SolutionSHA256'}),fullfile(outdir,'A_fair_test.csv'),'Encoding','UTF-8');
writetable(cell2table(summary,'VariableNames',{'Method','SelectedCandidate','Completed','Cases','SuccessRate','MeanItr', ...
    'MeanNF','MeanMedianTime','MaxFunctionNormSolved','MaxFeasibilitySolved','RepeatsVerified','StudyComplete'}), ...
    fullfile(outdir,'A_fair_summary.csv'),'Encoding','UTF-8');
writetable(cell2table(pairs,'VariableNames',{'Numerator','Denominator','CommonSuccesses','AllCases','NFGeometricRatio', ...
    'ItrGeometricRatio','TimeGeometricRatio','Scope'}),fullfile(outdir,'A_fair_paired.csv'),'Encoding','UTF-8');
end

function r=fair_georatio(a,b)
if isempty(a), r=NaN; else, r=exp(mean(log(a./b))); end
end

function fair_plots(fair,folder)
if ~isfolder(folder), mkdir(folder); end
labels=fair.protocol.labels;
exportOptions={'Units','centimeters','Width',18, ...
    'Padding','figure','PreserveAspectRatio','on'};
if all(fair.development.ConfigDone(1:5,:),'all')
    f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18 12.7]);
    cleanup=onCleanup(@()close(f)); %#ok<NASGU>
    boxes=[.065 .57 .27 .36;.38 .57 .27 .36;.695 .57 .27 .36; ...
        .2225 .105 .27 .36;.5375 .105 .27 .36];
    for c=1:5
        ax=axes(f,'Position',boxes(c,:));
        fair_profile(ax,squeeze(fair.development.NF(:,c,:)), ...
            squeeze(fair.development.Success(:,c,:)),labels);
        if c~=1 && c~=4, ylabel(ax,''); ax.YTickLabel=[]; end
        pbaspect(ax,[1.2 1 1]);
        title(ax,sprintf('(%c) \\sigma = %.4g','a'+c-1,fair.protocol.sigma_grid(c)), ...
            'FontWeight','normal','Interpreter','tex');
    end
    exportgraphics(f,fullfile(folder,'A_fair_sigma_profiles.pdf'),'ContentType','vector','Height',f.Position(4),exportOptions{:});
    exportgraphics(f,fullfile(folder,'A_fair_sigma_profiles.png'),'Resolution',300,'Height',f.Position(4),exportOptions{:});
    clear cleanup;
end
if fair.complete
    f=profile_figure();
    cleanup=onCleanup(@()close(f)); %#ok<NASGU>
    tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
    fields={'Itr','NF','Time'}; titles={'(a) Itr','(b) NF','(c) Tcpu'};
    for k=1:3
        ax=nexttile;
        fair_profile(ax,squeeze(fair.test.(fields{k})),squeeze(fair.test.Success),labels);
        if k~=1, ylabel(ax,''); end
        pbaspect(ax,[1.2 1 1]);
        title(ax,titles{k},'FontWeight','normal','FontSize',9);
    end
    exportgraphics(f,fullfile(folder,'A_fair_profiles.pdf'),'ContentType','vector','Height',f.Position(4),exportOptions{:});
    exportgraphics(f,fullfile(folder,'A_fair_profiles.png'),'Resolution',300,'Height',f.Position(4),exportOptions{:});
end
end

function fair_profile(ax,cost,success,labels)
colors=[0 114 178;0 158 115;204 121 167;213 94 0]/255;
styles={'-','-.',':','--'}; markers={'o','^','d','s'};
cost(~success | ~isfinite(cost) | cost<=0)=Inf;
best=min(cost,[],2); ratios=cost./best; ratios(~isfinite(ratios))=Inf;
finite=ratios(isfinite(ratios)); upper=2;
if ~isempty(finite), upper=2^ceil(log2(max(2,max(finite)*1.01))); end
profile_grid=unique([1;finite(:);upper]); hold(ax,'on');
for j=1:4
    y=arrayfun(@(q)mean(ratios(:,j)<=q*(1+1e-12)),profile_grid);
    stairs(ax,log2(profile_grid),y,'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',1.5);
    markerGrid=unique(round(linspace(1,numel(profile_grid),min(numel(profile_grid),20))));
    ix=markerGrid(j:4:end);
    plot(ax,log2(profile_grid(ix)),y(ix),markers{j},'Color',colors(j,:),'MarkerFaceColor','w', ...
        'MarkerSize',4,'LineStyle','none','HandleVisibility','off');
end
alphaUpper=log2(upper);
ticks=unique([0:max(1,ceil(alphaUpper/6)):alphaUpper,alphaUpper]);
set(ax,'XScale','linear','XTick',ticks,'FontName','Times New Roman','FontSize',8.5, ...
    'Box','on','TickDir','out','LineWidth',.7,'YTick',0:.25:1,'GridAlpha',.10, ...
    'XMinorGrid','off','YMinorGrid','off');
xlim(ax,[0 alphaUpper]); ylim(ax,[0 1.025]); grid(ax,'on');
xlabel(ax,'\alpha = log_2(\tau)'); ylabel(ax,'\rho_s(\alpha)');
legend(ax,labels,'Location','southeast','Box','off','Color','white','FontSize',7.2,'Interpreter','none');
end

function initial_step_study(options)
% Prespecified 2x2 ablation: warm initialization and one increase test.
o=struct('output_dir','','resume',true,'plots',true,'figs_only',false, ...
    'verify_only',false,'threads',1);
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'Unknown initial-step option %s',key{1});
    o.(key{1})=options.(key{1});
end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
if isempty(o.output_dir), o.output_dir=fullfile(root,'B组_步长策略'); end
outdir=fullfile(o.output_dir,'results'); figdir=fullfile(o.output_dir,'figures');
if ~isfolder(outdir), mkdir(outdir); end
if ~isfolder(figdir), mkdir(figdir); end
file=fullfile(outdir,'initial_steps.mat');
if o.figs_only
    d=load(file,'study'); initial_step_plot(d.study,figdir); return
end
validateattributes(o.threads,{'numeric'},{'scalar','integer','positive','finite'});
oldpath=path; oldrng=rng; oldthreads=maxNumCompThreads;
cleanup=onCleanup(@()abc_restore(oldpath,oldrng,oldthreads)); %#ok<NASGU>
maxNumCompThreads(o.threads);
src=source_dir(root,'monotone'); addpath(src,'-begin'); addpath(here,'-begin');
[~,p,meta]=load_cases(root,'monotone');
assert(abc_same_value(p,abc_main_parameters(meta.original_main,p.Itr_max)), ...
    'Stored TSITTCGA parameters differ from the current C driver.');
initial_step_verify(p);
if o.verify_only, fprintf('INITIAL_STEP_VERIFICATION_PASSED\n'); return; end
protocol=struct('date','2026-10-16','seed',20261016,'dimension',1000, ...
    'starts_per_problem',2,'problem_count',8,'initial_steps',[1e-4 .01 .1 .8], ...
    'warm_flags',[false false true true],'increase_flags',[false true false true], ...
    'timing_repeats',3,'max_iterations',2000,'epsilon',1e-6,'epsilon_C',1e-6, ...
    'threads',o.threads,'parameter_tuning',false, ...
    'note','Supplementary small-scale study, distinct from the formal 200 instances.');
labels={'Fixed','Fixed + increase','Warm start','Warm start + increase'};
q=struct('criterion','Proposed','sigma',p.sigma,'initial',p.rho_ls,'shrink',p.kappa, ...
    'warm',false,'increase',false,'lambda',.001,'nu',.8,'liu_mu',.01,'j_a',2,'j_l',.5,'label','');
configurations=cell(1,4);
for j=1:4
    c=q; c.warm=protocol.warm_flags(j); c.increase=protocol.increase_flags(j);
    c.label=labels{j}; configurations{j}=c;
end
controlled_source=fileread(fullfile(here,'TSITTCGA_LS.m'));
if o.resume && isfile(file)
    d=load(file,'study'); study=d.study;
    assert(abc_same_value(study.protocol,protocol) && abc_same_value(study.params,p) ...
        && strcmp(study.controlled_source,controlled_source) ...
        && strcmp(study.problem_source,meta.problem_source) ...
        && abc_same_value(study.configurations,configurations),'Initial-step protocol/source changed.');
    if all(study.Completed,'all')
        initial_step_tables(study,outdir);
        if o.plots, initial_step_plot(study,figdir); end
        fprintf('INITIAL_STEP_REUSE: all 256 configurations already completed.\n');
        return
    end
else
    rng(protocol.seed,'twister'); n=protocol.dimension;
    X0=cell(8,2);
    for id=1:8
        for r=1:2
            if id==3, X0{id,r}=2*rand(n,1)-1; else, X0{id,r}=4*rand(n,1); end
        end
    end
    cases=struct('Problem',{},'id',{},'Start',{},'n',{},'Delta',{});
    for delta=protocol.initial_steps
        for id=1:8
            for r=1:2
                cases(end+1)=struct('Problem',sprintf('P%d',id),'id',id, ...
                    'Start',r,'n',n,'Delta',delta); %#ok<AGROW>
            end
        end
    end
    M=nan(numel(cases),4);
    study=struct('protocol',protocol,'params',p,'configurations',{configurations}, ...
        'labels',{labels},'cases',cases,'X0',{X0}, ...
        'controlled_source',controlled_source,'problem_source',meta.problem_source, ...
        'driver_source',fileread([mfilename('fullpath'),'.m']),'matlab_version',version, ...
        'NF',M,'Itr',M,'LSNF',M,'Time',M,'RawNF',M,'RawTime',nan(numel(cases),4,3), ...
        'NormF',M,'Feasibility',M,'InitialAccepts',M,'Backtracks',M, ...
        'IncreaseTests',M,'IncreaseAccepts',M,'MedianStep',M, ...
        'Completed',false(size(M)),'Success',false(size(M)), ...
        'Status',{cell(size(M))},'SolutionHash',{cell(size(M))}, ...
        'AcceptedSteps',{cell(size(M))},'FirstProbePass',false(numel(cases),2));
    save(file,'study','-v7.3'); % Freeze the design and all starts before solving.
end
pf=str2func('problems');
assert(strcmpi(functions(pf).file,fullfile(src,'problems.m')),'Unexpected benchmark source.');
short=p; short.Itr_max=3;
pr=struct('name','P1','n',1000,'f',@(x,mode)pf(x,mode,1));
for j=1:4, TSITTCGA_LS(pr,study.X0{1,1},short,configurations{j}); end
clock=tic;
for i=1:numel(study.cases)
    c=study.cases(i); id=c.id; x0=study.X0{id,c.Start};
    pr=struct('name',c.Problem,'n',c.n,'f',@(x,mode)pf(x,mode,id));
    f0=pr.f(x0,1); d0=-f0;
    study.FirstProbePass(i,:)=[initial_step_test(x0,d0,c.Delta,pr,p), ...
        initial_step_test(x0,d0,min(1,c.Delta/p.kappa),pr,p)];
    % Probe evaluations above and independent final checks are outside solver NF/time.
    for j=mod(i-1+(0:3),4)+1
        if study.Completed(i,j), continue; end
        cfg=configurations{j}; cfg.initial=c.Delta;
        reference=[]; times=nan(1,3);
        for repeat=1:3
            [~,nf,it,~,info]=TSITTCGA_LS(pr,x0,p,cfg);
            times(repeat)=info.actual_time;
            normf=norm(pr.f(info.x,1)); feas=norm(info.x-pr.f(info.x,2));
            assert(info.actual_NF==1+info.inertia_nf+info.ls_nf+info.projection_nf, ...
                'Invalid total function evaluation accounting.');
            good=info.success && isfinite(normf) && normf<=1e-6 ...
                && isfinite(feas) && feas<=1e-6 && it<=p.Itr_max;
            assert(~info.success || good,'Invalid reported success.');
            assert(info.increase_tests<=info.initial_accepts ...
                && info.increase_accepts<=info.increase_tests,'Invalid increase-test counts.');
            if ~cfg.increase, assert(info.increase_tests==0); end
            if isempty(reference)
                reference=info; firstnf=nf; firstit=it; firstnorm=normf; firstfeas=feas;
            else
                assert(isequaln(info.x,reference.x) && info.actual_NF==reference.actual_NF ...
                    && isequaln(info.accepted_steps,reference.accepted_steps) ...
                    && info.success==reference.success,'Repeated solves changed numerical results.');
            end
        end
        info=reference; study.RawTime(i,j,:)=times;
        study.RawNF(i,j)=info.actual_NF; study.Success(i,j)=info.success;
        study.Status{i,j}=info.status; study.SolutionHash{i,j}=abc_hash(info.x);
        study.NormF(i,j)=firstnorm; study.Feasibility(i,j)=firstfeas;
        study.LSNF(i,j)=info.ls_nf; study.InitialAccepts(i,j)=info.initial_accepts;
        study.Backtracks(i,j)=info.backtracks; study.IncreaseTests(i,j)=info.increase_tests;
        study.IncreaseAccepts(i,j)=info.increase_accepts;
        study.AcceptedSteps{i,j}=info.accepted_steps;
        if ~isempty(info.accepted_steps), study.MedianStep(i,j)=median(info.accepted_steps); end
        if info.success
            study.NF(i,j)=firstnf; study.Itr(i,j)=firstit; study.Time(i,j)=median(times);
        end
        study.Completed(i,j)=true;
    end
    if mod(i,8)==0
        study.wall=toc(clock); save(file,'study','-v7.3');
        fprintf('INITIAL_STEP %d/%d cases, %.1fs\n',i,numel(study.cases),toc(clock));
    end
end
study.wall=toc(clock); save(file,'study','-v7.3');
initial_step_tables(study,outdir);
if o.plots, initial_step_plot(study,figdir); end
fprintf('INITIAL_STEP_DONE %d solves x %d timing repeats, %.1fs\n', ...
    numel(study.Success),protocol.timing_repeats,toc(clock));
end

function ok=initial_step_test(x,d,t,pr,p)
fz=pr.f(x+t*d,1); fw=pr.f(x,1);
factor=max(p.mu,norm(fz-fw)/(p.a+norm(fz)+norm(fw)));
ok=all(isfinite(fz)) && -fz'*d>=p.sigma*t*factor*norm(d)^2;
end

function initial_step_verify(p)
% Independent affine first-iteration check, with a manually computed projection.
pr=struct('name','AffineStepCheck','f',@initial_affine_problem);
q=p; q.Itr_max=1; x=ones(3,1); fw=2*x-1; direction=-fw;
cfg=struct('criterion','Proposed','sigma',p.sigma,'initial',1e-4,'shrink',p.kappa, ...
    'warm',false,'increase',false,'lambda',.001,'nu',.8,'liu_mu',.01,'j_a',2,'j_l',.5,'label','verify');
for j=1:2
    cfg.increase=j==2; t=cfg.initial;
    if cfg.increase, t=t/p.kappa; end
    assert(initial_step_test(x,direction,t,pr,p),'Analytic probe should pass.');
    z=x+t*direction; fz=2*z-1;
    expected=max(0,x-p.rho*(fz'*(x-z))/(fz'*fz)*fz);
    [~,~,~,~,info]=TSITTCGA_LS(pr,x,q,cfg);
    assert(norm(info.x-expected)<1e-13 && abs(info.accepted_steps(1)-t)<1e-15 ...
        && info.ls_nf==j && info.increase_tests==j-1 && info.increase_accepts==j-1 ...
        && info.initial_accepts==1 && info.backtracks==0,'Affine step/counter audit failed.');
end
fprintf('Independent affine update and extra-evaluation accounting passed.\n');
end

function value=initial_affine_problem(x,mode)
if mode==1, value=2*x-1; else, value=max(x,0); end
end

function initial_step_tables(s,folder)
rows=cell(0,19); summary=cell(0,14);
for delta=s.protocol.initial_steps
    mask=[s.cases.Delta]'==delta;
    for j=1:4
        ok=s.Success(mask,j); tests=sum(s.IncreaseTests(mask,j));
        accept_rate=NaN; if tests>0, accept_rate=sum(s.IncreaseAccepts(mask,j))/tests; end
        summary(end+1,:)={delta,s.labels{j},sum(mask),mean(ok),mean(s.Itr(mask,j)), ...
            mean(s.NF(mask,j)),complete_mean(s.LSNF(mask,j),ok),mean(s.Time(mask,j)), ...
            mean(s.IncreaseTests(mask,j)),accept_rate,mean(s.FirstProbePass(mask,1)), ...
            mean(all(s.FirstProbePass(mask,:),2)),mean(s.RawNF(mask,j)), ...
            mean(s.NF(mask & s.Success(:,j),j))}; %#ok<AGROW>
        for i=find(mask)'
            c=s.cases(i);
            rows(end+1,:)={c.Problem,c.n,c.Start,delta,s.labels{j},s.Success(i,j), ...
                s.Itr(i,j),s.NF(i,j),s.LSNF(i,j),s.Time(i,j),s.NormF(i,j), ...
                s.Feasibility(i,j),s.InitialAccepts(i,j),s.Backtracks(i,j), ...
                s.IncreaseTests(i,j),s.IncreaseAccepts(i,j),s.FirstProbePass(i,1), ...
                s.FirstProbePass(i,2),s.Status{i,j}}; %#ok<AGROW>
        end
    end
end
T=cell2table(rows,'VariableNames',{'Problem','n','Start','InitialStep','Strategy', ...
    'Success','Itr','NF','LSNF','Time','NormF','Feasibility','InitialAccepts','Backtracks', ...
    'IncreaseTests','IncreaseAccepts','FirstStepPass','FirstEnlargedStepPass','Status'});
writetable(T,fullfile(folder,'initial_steps_runs.csv'),'Encoding','UTF-8');
T=cell2table(summary,'VariableNames',{'InitialStep','Strategy','Instances','SuccessRate', ...
    'MeanItr','MeanNF','MeanLSNF','MeanTime','MeanIncreaseTests','IncreaseAcceptanceRate', ...
    'FirstStepPassRate','BothFirstProbesPassRate','MeanRawNF','MeanNFSolved'});
writetable(T,fullfile(folder,'initial_steps_summary.csv'),'Encoding','UTF-8');
disp(T(:,[1 2 4 6 7 10]));
end

function initial_step_plot(s,folder)
colors=[70 70 70;213 94 0;0 158 115;0 114 178]/255;
styles={'--',':','-.','-'}; markers={'s','d','^','o'};
f=profile_figure();
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
for m=1:3
    ax=nexttile; hold(ax,'on');
    if m==1
        for j=1:4
            y=arrayfun(@(delta)mean(s.Success([s.cases.Delta]'==delta,j)),s.protocol.initial_steps);
            plot(ax,s.protocol.initial_steps,y,'Color',colors(j,:),'LineStyle',styles{j}, ...
                'Marker',markers{j},'LineWidth',1.4,'MarkerFaceColor','w','MarkerSize',4);
        end
        set(ax,'XScale','log','XTick',s.protocol.initial_steps);
        xlim(ax,[s.protocol.initial_steps(1)*.7 s.protocol.initial_steps(end)*1.4]);
        xlabel(ax,'Initial trial step'); ylabel(ax,'Fraction solved');
        title(ax,'(a) Initial-step sensitivity','FontWeight','normal','FontSize',9);
    else
        cost=s.NF; label='Function evaluations';
        if m==3, cost=s.LSNF; label='Line search evaluations'; end
        cost(~s.Success)=Inf; ratios=cost./min(cost,[],2); ratios(~isfinite(ratios))=Inf;
        finite=ratios(isfinite(ratios)); upper=2;
        if ~isempty(finite), upper=2^ceil(log2(max(2,max(finite)*1.01))); end
        gridx=unique([1;finite(:);upper]);
        for j=1:4
            y=arrayfun(@(v)mean(ratios(:,j)<=v*(1+1e-12)),gridx);
            stairs(ax,log2(gridx),y,'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',1.5);
            ix=unique(round(linspace(1,numel(gridx),min(20,numel(gridx))))); ix=ix(j:4:end);
            plot(ax,log2(gridx(ix)),y(ix),markers{j},'Color',colors(j,:), ...
                'MarkerSize',4,'LineStyle','none','MarkerFaceColor','w','HandleVisibility','off');
        end
        xmax=log2(upper); xlim(ax,[0 xmax]); set(ax,'XTick',unique([0:max(1,ceil(xmax/6)):xmax,xmax]));
        xlabel(ax,'\alpha = log_2(\tau)');
        title(ax,sprintf('(%c) %s','a'+m-1,label),'FontWeight','normal','FontSize',9);
    end
    set(ax,'FontName','Times New Roman','FontSize',8.5,'Box','on','TickDir','out', ...
        'LineWidth',.7,'YTick',0:.25:1,'GridAlpha',.10,'XMinorGrid','off','YMinorGrid','off');
    ylim(ax,[0 1.025]); grid(ax,'on'); pbaspect(ax,[1.2 1 1]);
    legend(ax,s.labels,'Location','southeast','Box','off','FontSize',6.5,'Interpreter','none');
end
args={'Units','centimeters','Width',f.Position(3),'Height',f.Position(4),'Padding','figure','PreserveAspectRatio','on'};
exportgraphics(f,fullfile(folder,'initial_steps_profiles.pdf'),'ContentType','vector',args{:});
exportgraphics(f,fullfile(folder,'initial_steps_profiles.png'),'Resolution',300,args{:});
end

function step_probe_study(options)
% Matched first-search diagnostics, executed by the unchanged controlled solver.
o=struct('output_dir','','resume',true,'plots',true,'figs_only',false, ...
    'verify_only',false,'threads',1);
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'Unknown step-probe option %s',key{1});
    o.(key{1})=options.(key{1});
end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
if isempty(o.output_dir), o.output_dir=fullfile(root,'B组_步长策略'); end
outdir=fullfile(o.output_dir,'results'); figdir=fullfile(o.output_dir,'figures');
file=fullfile(outdir,'initial_steps.mat');
assert(isfile(file),'Run the existing initial_steps study first.');
saved=load(file); study=saved.study;
if o.figs_only
    assert(isfield(saved,'local_probe'),'Run mode=step_probe first.');
    step_probe_outputs(saved.local_probe,outdir,figdir,true); return
end
validateattributes(o.threads,{'numeric'},{'scalar','integer','positive','finite'});
oldpath=path; oldrng=rng; oldthreads=maxNumCompThreads;
cleanup=onCleanup(@()abc_restore(oldpath,oldrng,oldthreads)); %#ok<NASGU>
maxNumCompThreads(o.threads);
src=source_dir(root,'monotone'); addpath(src,'-begin'); addpath(here,'-begin');
clear problems; pf=str2func('problems');
assert(strcmpi(functions(pf).file,fullfile(src,'problems.m')),'Wrong benchmark source.');
assert(strcmpi(which('TSITTCGA_LS'),fullfile(here,'TSITTCGA_LS.m')),'Wrong controlled solver.');
p=study.params;
assert(abc_same_value(p,abc_main_parameters(fileread(fullfile(src,'main.m')),p.Itr_max)), ...
    'Stored parameters differ from the current C driver.');
assert(strcmp(study.controlled_source,fileread(fullfile(here,'TSITTCGA_LS.m'))) ...
    && strcmp(study.problem_source,fileread(fullfile(src,'problems.m'))), ...
    'Stored experiment sources have changed.');
initial_step_verify(p);
if o.verify_only, fprintf('STEP_PROBE_VERIFICATION_PASSED\n'); return; end
protocol=struct('version',1,'date','2026-10-16','seed',study.protocol.seed, ...
    'dimension',study.protocol.dimension,'starts_per_problem',2,'problem_count',8, ...
    'initial_steps',unique([logspace(-6,0,121),study.protocol.initial_steps]), ...
    'first_search_only',true,'xbar_equals_x0',true,'warm_start',false, ...
    'increase_flags',[false true],'parameter_tuning',false, ...
    'note','Paired local search outcomes, not whole-solve performance or independent samples.');
if o.resume && isfield(saved,'local_probe')
    local_probe=saved.local_probe;
    assert(abc_same_value(local_probe.protocol,protocol) && abc_same_value(local_probe.params,p) ...
        && isequaln(local_probe.X0,study.X0),'Stored local-probe design differs.');
else
    cfg=study.configurations{1}; cfg.warm=false;
    q=p; q.Itr_max=1; q.quiet=true;
    rows=cell(numel(protocol.initial_steps)*16,14); row=0;
    clock=tic;
    for id=1:8
        pr=struct('name',sprintf('P%d',id),'n',protocol.dimension,'f',@(x,mode)pf(x,mode,id));
        for r=1:2
            x0=study.X0{id,r}; q.xbar=x0; fw=pr.f(x0,1); direction=-fw;
            assert(norm(fw)>1e-6,'A local probe must start before convergence.');
            for delta=protocol.initial_steps
                cfg.initial=delta; pair=cell(1,2);
                for j=1:2
                    cfg.increase=j==2;
                    [~,~,~,~,v]=TSITTCGA_LS(pr,x0,q,cfg);
                    assert(numel(v.accepted_steps)==1,'The first search did not complete.');
                    assert(v.actual_NF==1+v.inertia_nf+v.ls_nf+v.projection_nf ...
                        && v.ls_nf==1+v.backtracks+v.increase_tests,'Incorrect NF accounting.');
                    assert(initial_step_test(x0,direction,v.accepted_steps(1),pr,p), ...
                        'Returned step fails the acceptance criterion.');
                    pair{j}=v;
                end
                off=pair{1}; on=pair{2}; t0=off.accepted_steps(1); t1=on.accepted_steps(1);
                assert(off.initial_accepts==on.initial_accepts && off.backtracks==on.backtracks);
                assert(off.increase_tests==0 && on.ls_nf-off.ls_nf==on.increase_tests);
                if off.initial_accepts==0
                    outcome='Initial rejected';
                    assert(on.increase_tests==0 && t0==t1 && off.backtracks>=1);
                    assert(~initial_step_test(x0,direction,t0/p.kappa,pr,p), ...
                        'The preceding backtracking trial should have failed.');
                elseif on.increase_tests==0
                    outcome='No increase needed';
                    assert(delta==1 || norm(pr.f(x0+delta*direction,1))==0);
                elseif on.increase_accepts==1
                    outcome='Both accepted';
                    assert(t0==delta && t1==min(1,delta/p.kappa) && t1>t0);
                else
                    outcome='Larger rejected';
                    assert(t0==delta && t1==t0 ...
                        && ~initial_step_test(x0,direction,min(1,delta/p.kappa),pr,p));
                end
                if t0==t1, assert(isequaln(off.x,on.x),'Same step must give the same update.'); end
                row=row+1;
                rows(row,:)={id,r,protocol.dimension,delta,outcome,off.initial_accepts, ...
                    on.increase_tests,on.increase_accepts,off.backtracks, ...
                    t0,t1,off.ls_nf,on.ls_nf,t1/t0};
            end
        end
        fprintf('STEP_PROBE P%d/8 completed\n',id);
    end
    T=cell2table(rows,'VariableNames',{'Problem','Start','n','InitialStep','Outcome', ...
        'InitialAccepted','IncreaseTests','IncreaseAccepted','Backtracks', ...
        'StepWithoutIncrease','StepWithIncrease','LSNFWithoutIncrease','LSNFWithIncrease','StepRatio'});
    local_probe=struct('protocol',protocol,'params',p,'X0',{study.X0},'table',T, ...
        'controlled_source',study.controlled_source,'problem_source',study.problem_source, ...
        'driver_source',fileread([mfilename('fullpath'),'.m']),'elapsed',toc(clock));
    save(file,'local_probe','-append'); % Preserve the pre-existing global study variable.
    check=load(file,'study'); assert(abc_same_value(check.study,study),'The global study changed.');
end
step_probe_outputs(local_probe,outdir,figdir,o.plots);
fprintf('STEP_PROBE_DONE: %d paired searches on 16 fixed initial states.\n',height(local_probe.table));
end

function step_probe_outputs(s,outdir,figdir,plots)
T=s.table; writetable(T,fullfile(outdir,'initial_steps_probe.csv'),'Encoding','UTF-8');
labels={'Both accepted','Larger rejected','Initial rejected','No increase needed'};
for j=1:4
    mask=strcmp(T.Outcome,labels{j});
    fprintf('%s: %d pairs, mean extra LSNF %.3g, mean step ratio %.3g\n', ...
        labels{j},sum(mask),mean(T.LSNFWithIncrease(mask)-T.LSNFWithoutIncrease(mask)),mean(T.StepRatio(mask)));
end
if ~plots, return; end
if ~isfolder(figdir), mkdir(figdir); end
gridx=s.protocol.initial_steps; share=zeros(numel(gridx),4); steps=zeros(numel(gridx),2); nf=steps;
for i=1:numel(gridx)
    mask=T.InitialStep==gridx(i); assert(sum(mask)==16);
    for j=1:4, share(i,j)=mean(strcmp(T.Outcome(mask),labels{j})); end
    steps(i,:)=[mean(T.StepWithoutIncrease(mask)),mean(T.StepWithIncrease(mask))];
    nf(i,:)=[mean(T.LSNFWithoutIncrease(mask)),mean(T.LSNFWithIncrease(mask))];
end
f=profile_figure();
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
colors=[0 114 178;213 94 0;70 70 70;0 158 115]/255;
styles={'-','--',':','-.'};
for m=1:3
    ax=nexttile; hold(ax,'on');
    if m==1
        active=find(any(share>0,1));
        for j=active
            stairs(ax,gridx,share(:,j),'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',1.3);
        end
        ylabel(ax,'Fraction of initial states'); ylim(ax,[0 1.03]);
        title(ax,'(a) First-trial outcomes','FontWeight','normal');
        legend(ax,labels(active),'Location','west','FontSize',6.7,'Box','off');
    else
        values=steps; ylab='Mean accepted step'; cap='(b) Accepted steps';
        if m==3, values=nf; ylab='Mean line search evaluations'; cap='(c) Evaluation cost'; end
        plot(ax,gridx,values(:,1),'--','Color',colors(3,:),'LineWidth',1.5);
        plot(ax,gridx,values(:,2),'-','Color',colors(1,:),'LineWidth',1.5);
        ylabel(ax,ylab); title(ax,cap,'FontWeight','normal');
        if m==2, set(ax,'YScale','log'); else, ylim(ax,[0 max(nf,[],'all')*1.1]); end
        legend(ax,{'Without increase','With increase'},'Location','northwest','FontSize',7,'Box','off');
    end
    set(ax,'XScale','log','XTick',[1e-6 1e-4 1e-2 1],'XLim',[1e-6 1], ...
        'FontName','Times New Roman','FontSize',8.5,'Box','on','TickDir','out', ...
        'LineWidth',.7,'GridAlpha',.10,'XMinorGrid','off','YMinorGrid','off');
    xlabel(ax,'Initial trial step'); grid(ax,'on'); pbaspect(ax,[1.2 1 1]);
end
args={'Units','centimeters','Width',f.Position(3),'Height',f.Position(4),'Padding','figure','PreserveAspectRatio','on'};
exportgraphics(f,fullfile(figdir,'initial_steps_probe.pdf'),'ContentType','vector',args{:});
exportgraphics(f,fullfile(figdir,'initial_steps_probe.png'),'Resolution',300,args{:});
end













