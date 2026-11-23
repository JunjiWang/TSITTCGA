function T=run_pmu(mode,systems,budget,repetitions)
% All systems use the projector backend fixed by pmu_parameters.m.
%RUN_PMU Local, standalone comparison on frozen IEEE14/118/300 PMU data.
% With no argument, run all selected cases in final mode.
% run_pmu('check') : all 45 IEEE14 cases, cases 1/16/31 of IEEE118 and300, 1 repeat.
% run_pmu('final') : IEEE14 + IEEE118, 3 timed repeats per solver.
% Optional arguments override systems, common iteration budget, and repeats.
here=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;
oldthreads=maxNumCompThreads(1);
cleanup=onCleanup(@()restore_session(oldpath,oldrng,oldthreads)); %#ok<NASGU>
addpath(here,'-begin');
helpers={'pmu_parameters','pmu_map','pmu_run_cgpm','pmu_prepare_backend', ...
    'huber_irls_same_model','summarize_pmu'};
for h=helpers
    assert(strcmpi(fileparts(which(h{1})),here),'PMU:Shadowing', ...
        'A different %s is active. Change the current folder to this package.',h{1});
end
c=pmu_parameters();
if nargin<1,mode='final';end
assert(any(strcmp(mode,{'check','final'})),'Use check or final.');
if nargin<2||isempty(systems)
    if strcmp(mode,'check'),systems=c.systems;else,systems=c.main_systems;end
end
if ischar(systems)||isstring(systems),systems=cellstr(systems);end
if nargin<3||isempty(budget),budget=c.budget;end
if nargin<4||isempty(repetitions)
    repetitions=c.repetitions;if strcmp(mode,'check'),repetitions=1;end
end
assert(budget>=1&&budget==fix(budget)&&repetitions>=1&&repetitions==fix(repetitions));
assert(all(ismember(systems,c.systems)));
for a=1:numel(c.solvers)
    f=functions(c.solvers{a});assert(strcmpi(fileparts(f.file),here),'Wrong solver source: %s',f.file);
    c.params{a}.Itr_max=budget;
end
out=fullfile(here,'results');if ~isfolder(out),mkdir(out);end
prefix=sprintf('%s_budget%d_rep%d',mode,budget,repetitions);
if ~isequal(systems,c.systems),prefix=[prefix '_' strjoin(systems,'_')];end
prefix=[prefix '_' c.profile];
matfile=fullfile(out,[prefix '.mat']);
assert(~isfile(matfile),'PMU:ExistingResults',['Results already exist: %s. ' ...
        'Preserve them before rerunning; main archives old results for a full rerun.'],matfile);
protocol=struct('Created',char(datetime('now')),'MATLAB',version,'Mode',mode, ...
    'Budget',budget,'Repetitions',repetitions,'Config',c,'Threads',1, ...
    'Timing','elapsed tic/toc, common preparation separate, solver and total both saved', ...
    'Initialization','CGPM: positive/negative split of P*y; IRLS: OLS state', ...
    'NoModelScaling',true,'NoRestarts',true,'NoEndpointProjection',true);
protocol.SelectedSystems=systems;
protocol.ParameterProfile=c.profile;
protocol.InitialPoint=c.initial;
protocol.FullCoverage=strcmp(mode,'final');
files=dir(fullfile(here,'*.m'));snapshot=struct([]);
for k=1:numel(files)
    snapshot(k).Name=files(k).name; %#ok<AGROW>
    snapshot(k).Text=fileread(fullfile(here,files(k).name)); %#ok<AGROW>
end
protocol.SourceSnapshot=snapshot;
rows=struct([]);endpoints={};preparation=struct([]);completed=false;
for sy=1:numel(systems)
    sys=systems{sy};d=load(fullfile(here,'data',['input_' sys '.mat']));
    if isfield(d,'instance_records'),cases=d.instance_records;else,cases=d.instances;end
    if strcmp(mode,'check')&&~strcmp(sys,'ieee14'),idx=[1 16 31];else,idx=1:numel(cases);end
    % JIT warm-up, never included in timed samples or endpoint success rates.
    backend=instance_backend(c,sys);
[warm,~]=pmu_prepare_backend(cases{idx(1)},backend);
    for a=1:6
        if a<=5
            p=c.params{a};p.Itr_max=3;pmu_run_cgpm(warm,c.solvers{a},p,c.lambda,c.initial);
        else,huber_irls_same_model(warm,c.lambda,c.target,3,false);end
    end
    for j=1:numel(idx)
        k=idx(j);[inst,prep]=pmu_prepare_backend(cases{k},backend);
        prep.System=sys;prep.Case=k;
        if isempty(preparation),preparation=prep;else,preparation(end+1)=prep;end %#ok<AGROW>
        for rep=1:repetitions
            order=mod((0:5)+(j-1)+(rep-1),6)+1;
            for ai=1:6
                a=order(ai);
                if a<=5
                    r=pmu_run_cgpm(inst,c.solvers{a},c.params{a},c.lambda,c.initial);
                else
                    r=huber_irls_same_model(inst,c.lambda,c.target,budget,false);
                    r.ActualItr=r.Itr;r.ActualNF=NaN;r.count_ok=true;r.report_ok=true;
                    r.PaperFloorAccepts=0;r.criterion_reached=r.residual<=c.target&&r.feasibility<=c.feasibility_target;
                    % Validate IRLS using the same independent state-LS identity.
                    fit=inst.y-r.e-inst.G*r.xi;
                    fr=[min(r.s(1:inst.m),c.lambda-fit);min(r.s(inst.m+1:end),c.lambda+fit)];
                    assert(abs(norm(fr)-r.residual)<=1e-9,'IRLS independent residual mismatch.');
                    assert(r.monotonicity_violations==0,'IRLS objective increased.');
                end
                row=make_row(sys,k,rep,c.names{a},inst,r);
                if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
                endpoints{end+1}=r; %#ok<AGROW>
                fprintf('%s case %2d rep %d %-14s ok=%d itr=%5d F=%.3e total=%.3fs\n', ...
                    sys,k,rep,c.names{a},r.success,r.ActualItr,r.residual,r.TotalTime);
            end
        end
        save(matfile,'rows','endpoints','preparation','protocol','completed','-v7.3');
        writetable(struct2table(rows),fullfile(out,[prefix '_runs.csv']));
    end
end
completed=true;save(matfile,'rows','endpoints','preparation','protocol','completed','-v7.3');
T=summarize_pmu(matfile);disp(T);
end
function row=make_row(sys,k,rep,name,inst,r)
row=struct('System',sys,'Case',k,'Repeat',rep,'Method',name,'Success',r.success, ...
    'CriterionReached',r.criterion_reached,'Status',r.status,'Measurements',inst.m, ...
    'StateDimension',size(inst.G,2),'CNEDimension',2*inst.m,'BadData',nnz(inst.bad), ...
    'Itr',r.ActualItr,'NF',r.ActualNF,'Residual',r.residual,'Feasibility',r.feasibility, ...
    'Objective',r.obj,'RelativeStateError',r.relerr,'SolverWallTime',r.SolverWallTime, ...
    'TotalTime',r.TotalTime,'CountMatch',r.count_ok,'ReportMatch',r.report_ok, ...
    'WeightedLSSolves',r.WeightedLSSolves,'CNEChecks',r.CNEChecks,'PaperFloorAccepts',r.PaperFloorAccepts);
end
function restore_session(p,r,t)
path(p);rng(r);maxNumCompThreads(t);
end



function b=instance_backend(c,sys)
% Per-system projector backend: the two implementations compute the same
% orthogonal projector and differ only by machine-precision rounding.
if isfield(c,'backend_by_system') && isfield(c.backend_by_system,sys)
    b=c.backend_by_system.(sys);
else
    b=c.backend;
end
end
