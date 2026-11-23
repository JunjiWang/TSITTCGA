function [T,S]=run_pmu_qp_validation()
%RUN_PMU_QP_VALIDATION Independent accuracy checks on ten fixed PMU cases.
%   No parameter tuning, replacement model, endpoint repair, or speed ranking.
%   quadprog solves the sparse lifted QP at its stricter native tolerances;
%   both methods are independently checked against the original 1e-6 CNE.
%   The selection below is a fixed coverage sample, not a random sample and
%   not a claim that the optimum or the estimated state is unique.
%
%   Files: results/qp_validation_<parameter profile>.mat,
%          ..._selection.csv, ..._cases.csv, and ..._summary.csv.
%   Existing output is never overwritten. Checkpoints retain failed solves.
here=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;
oldthreads=maxNumCompThreads(1);
cleanup=onCleanup(@()restore_session(oldpath,oldrng,oldthreads)); %#ok<NASGU>
addpath(here,'-begin');
required={'pmu_parameters','pmu_prepare_backend','pmu_run_cgpm','pmu_map', ...
    'pmu_qp_reference','TSITTCGA_precision'};
for k=1:numel(required)
    assert(strcmpi(fileparts(which(required{k})),here),'PMUQP:Shadowing', ...
        'A different %s shadows the source in this package.',required{k});
end
assert(exist('quadprog','file')~=0&&exist('optimoptions','file')~=0 ...
    &&license('test','Optimization_Toolbox'),'PMUQP:MissingToolbox', ...
    'The independent QP reference requires MATLAB Optimization Toolbox.');
c=pmu_parameters();
% Stricter reference tolerances are applied uniformly, never selected by the
% TS result. The initial native-1e-10 study remains in separate result files.
reference_options=struct('NativeTolerance',1e-14,'StepTolerance',1e-16);
systems={'ieee14','ieee118','ieee300'};
case_ids={[1 23 45],[3 22 41],[5 23 31 45]};
out=fullfile(here,'results');if ~isfolder(out),mkdir(out);end
prefix=['qp_validation_' c.profile];matfile=fullfile(out,[prefix '.mat']);
suffixes={'.mat','_selection.csv','_cases.csv','_summary.csv'};
for k=1:numel(suffixes)
    assert(~isfile(fullfile(out,[prefix suffixes{k}])),'PMUQP:ExistingResults', ...
        'Output already exists: %s. Preserve it before rerunning.', ...
        fullfile(out,[prefix suffixes{k}]));
end

% Freeze the selection, input identities and source before observing results.
selection_rows=struct([]);input_files=struct([]);
for sy=1:numel(systems)
    input_path=fullfile(here,'data',['input_' systems{sy} '.mat']);
    d=load(input_path,'instances');cases=d.instances;
    input_files(sy).Name=['input_' systems{sy} '.mat']; %#ok<AGROW>
    input_files(sy).SHA256=sha256_file(input_path); %#ok<AGROW>
    for k=case_ids{sy}
        a=cases{k};
        row=struct('SelectionIndex',numel(selection_rows)+1,'System',systems{sy}, ...
            'Case',k,'Measurements',a.m,'StateDimension',size(a.G,2), ...
            'CNEDimension',2*a.m,'Currents',a.Currents, ...
            'BadData',nnz(a.bad),'ConfiguredBadData',a.BadData, ...
            'SeedIndex',a.SeedIndex,'RandomSeed',a.RandomSeed);
        selection_rows=push(selection_rows,row);
    end
end
selection=struct2table(selection_rows);
assert(height(selection)==10,'PMUQP:Selection','The fixed selection must contain ten cases.');
files=dir(fullfile(here,'*.m'));snapshot=struct([]);
for k=1:numel(files)
    p=fullfile(here,files(k).name);
    snapshot(k).Name=files(k).name;snapshot(k).Text=fileread(p); %#ok<AGROW>
    snapshot(k).SHA256=sha256_file(p); %#ok<AGROW>
end
protocol=struct('Created',char(datetime('now','TimeZone','Asia/Shanghai')), ...
    'MATLAB',version,'Purpose','Independent convex-QP numerical accuracy reference', ...
    'SelectionRule','Fixed in source before solves; coverage sample, not random or result-selected', ...
    'SelectedSystems',{systems},'SelectedCaseIDs',{case_ids},'Config',c, ...
    'ParameterProfile',c.profile,'Threads',1,'OriginalTarget',c.target, ...
    'FeasibilityTarget',c.feasibility_target,'TSBudget',c.budget, ...
    'TSInitialization',c.initial,'QPNativeTolerance',reference_options.NativeTolerance, ...
    'QPNativeStepTolerance',reference_options.StepTolerance,'ReferenceOptions',reference_options, ...
    'QPReferenceOnly',true,'SamePrecisionTimingComparison',false, ...
    'TimingInterpretation','Single-run reference metadata only; no efficiency ranking or median-time claim', ...
    'UniqueOptimumClaim',false,'ExactCertificate',false,'NoTuning',true, ...
    'NoRestarts',true,'NoEndpointRepair',true,'FailureRecordsRetained',true, ...
    'InputFiles',input_files,'SourceSnapshot',snapshot);
rows=struct([]);endpoints=cell(10,2);fresh_checks=cell(10,2);
preparation=cell(10,1);completed=false;
save(matfile,'protocol','selection','rows','endpoints','fresh_checks', ...
    'preparation','completed','-v7.3');
writetable(selection,fullfile(out,[prefix '_selection.csv']));

for j=1:height(selection)
    sys=selection.System{j};k=selection.Case(j);
    d=load(fullfile(here,'data',['input_' sys '.mat']),'instances');
    backend=c.backend;
    if isfield(c,'backend_by_system'),backend=c.backend_by_system.(sys);end
    [inst,prep]=pmu_prepare_backend(d.instances{k},backend);
    preparation{j}=prep;
    % Neither solution is used as the other method's initialization.
    q=pmu_qp_reference(inst,c.lambda,c.target,reference_options);
    a=pmu_run_cgpm(inst,c.solvers{1},c.params{1},c.lambda,c.initial);
    qfresh=fresh_check(q,inst,c);afresh=fresh_check(a,inst,c);
    endpoints{j,1}=a;endpoints{j,2}=q;
    fresh_checks{j,1}=afresh;fresh_checks{j,2}=qfresh;
    % An inconsistent finite report is an implementation error. A solver
    % failure or an unmet target remains an ordinary, retained outcome.
    validate_report(a,afresh,'TSITTCGA');validate_report(q,qfresh,'quadprog');
    ts_verified=logical(a.success)&&afresh.TargetReached&&a.count_ok&&a.report_ok;
    qp_verified=logical(q.success)&&qfresh.TargetReached;
    absdiff=abs(afresh.Objective-qfresh.Objective);
    reldiff=absdiff/max(abs(qfresh.Objective),realmin);
    scaled_diff=absdiff/max([1,abs(afresh.Objective),abs(qfresh.Objective)]);
    state_distance=NaN;
    if afresh.Finite&&qfresh.Finite
        state_distance=norm(afresh.State-qfresh.State)/max(norm(qfresh.State),realmin);
    end
    row=struct('System',sys,'Case',k,'Measurements',inst.m, ...
        'StateDimension',size(inst.G,2),'CNEDimension',2*inst.m, ...
        'Currents',selection.Currents(j),'BadData',nnz(inst.bad), ...
        'RandomSeed',selection.RandomSeed(j),'TSVerified',ts_verified, ...
        'QPVerified',qp_verified,'TSStatus',a.status,'QPStatus',q.status, ...
        'TSOriginalResidual',afresh.Residual,'QPOriginalResidual',qfresh.Residual, ...
        'TSFeasibility',afresh.Feasibility,'QPFeasibility',qfresh.Feasibility, ...
        'TSObjective',afresh.Objective,'QPObjective',qfresh.Objective, ...
        'SignedObjectiveDifference',afresh.Objective-qfresh.Objective, ...
        'AbsoluteObjectiveDifference',absdiff,'RelativeObjectiveDifference',reldiff, ...
        'ScaledObjectiveDifference',scaled_diff, ...
        'TSRelativeStateError',afresh.RelativeStateError, ...
        'QPRelativeStateError',qfresh.RelativeStateError, ...
        'RelativeStateVectorDifference',state_distance,'TSItr',a.ActualItr, ...
        'TSNF',a.ActualNF,'TSNFCounted',a.NFCounted,'TSCountMatch',a.count_ok, ...
        'QPExitflag',q.exitflag,'QPItr',q.Itr, ...
        'QPKKTScaled',nested_value(q,'KKT','ResidualScaled'), ...
        'QPDualityGap',nested_value(q,'Duality','Gap'), ...
        'QPScaledDualityGap',nested_value(q,'Duality','ScaledGap'), ...
        'TSFreshResidualDifference',afresh.ReportedResidualDifference, ...
        'QPFreshResidualDifference',qfresh.ReportedResidualDifference, ...
        'TSFreshObjectiveDifference',afresh.ReportedObjectiveDifference, ...
        'QPFreshObjectiveDifference',qfresh.ReportedObjectiveDifference, ...
        'TSReferenceRunSolverWallTime',a.SolverWallTime, ...
        'QPReferenceRunSolverWallTime',q.SolverWallTime, ...
        'TSReferenceRunTotalTime',a.TotalTime,'QPReferenceRunTotalTime',q.TotalTime);
    rows=push(rows,row);
    fprintf('QP validation %2d/10 %s case %2d: TS=%d QP=%d F=(%.3e,%.3e) |df|=%.3e rel=%.3e\n', ...
        j,sys,k,ts_verified,qp_verified,afresh.Residual,qfresh.Residual,absdiff,reldiff);
    save(matfile,'protocol','selection','rows','endpoints','fresh_checks', ...
        'preparation','completed','-v7.3');
    writetable(struct2table(rows),fullfile(out,[prefix '_cases.csv']));
end
T=struct2table(rows);S=summarize(T,systems);completed=true;
save(matfile,'protocol','selection','rows','endpoints','fresh_checks', ...
    'preparation','completed','T','S','-v7.3');
writetable(S,fullfile(out,[prefix '_summary.csv']));
disp(S);
end

function a=fresh_check(r,inst,c)
% Fresh sparse backslash does not reuse the prepared QR recovery object.
a=struct('Finite',false,'TargetReached',false,'State',[], ...
    'Residual',Inf,'Feasibility',Inf,'Objective',NaN,'RelativeStateError',NaN, ...
    'ReportedResidualDifference',NaN,'ReportedObjectiveDifference',NaN);
if ~isfield(r,'s')||numel(r.s)~=2*inst.m||any(~isfinite(r.s)),return;end
s=r.s(:);e=s(1:inst.m)-s(inst.m+1:end);G=sparse(inst.G);y=inst.y(:);
xi=G\(y-e);fit=y-e-G*xi;
f=[min(s(1:inst.m),c.lambda-fit);min(s(inst.m+1:end),c.lambda+fit)];
a.State=xi;a.Residual=norm(f);a.Feasibility=norm(min(s,0));
a.Objective=.5*dot(fit,fit)+c.lambda*norm(e,1);
a.RelativeStateError=norm(xi-inst.true_xi)/norm(inst.true_xi);
a.Finite=all(isfinite([xi;f;a.Objective;a.RelativeStateError]));
a.TargetReached=a.Finite&&a.Residual<=c.target&&a.Feasibility<=c.feasibility_target;
if isfinite(r.residual),a.ReportedResidualDifference=abs(a.Residual-r.residual);end
if isfinite(r.obj),a.ReportedObjectiveDifference=abs(a.Objective-r.obj);end
end

function validate_report(r,a,name)
if ~a.Finite,return;end
if isfinite(r.residual)
    tol=1e-9*max([1,abs(r.residual),a.Residual]);
    assert(a.ReportedResidualDifference<=tol,'PMUQP:ResidualMismatch', ...
        '%s residual disagrees with fresh sparse-LS recomputation.',name);
end
if isfinite(r.obj)
    tol=1e-10*max([1,abs(r.obj),abs(a.Objective)]);
    assert(a.ReportedObjectiveDifference<=tol,'PMUQP:ObjectiveMismatch', ...
        '%s objective disagrees with fresh sparse-LS recomputation.',name);
end
end

function S=summarize(T,systems)
rows=struct([]);groups=[systems,{'all'}];
for k=1:numel(groups)
    if strcmp(groups{k},'all'),ix=true(height(T),1);else,ix=strcmp(T.System,groups{k});end
    z=T(ix,:);common=z.TSVerified&z.QPVerified;
    row=struct('System',groups{k},'Cases',height(z),'TSVerified',nnz(z.TSVerified), ...
        'QPVerified',nnz(z.QPVerified),'CommonVerified',nnz(common), ...
        'MaxTSResidualVerified',finite_max(z.TSOriginalResidual(z.TSVerified)), ...
        'MaxQPResidualVerified',finite_max(z.QPOriginalResidual(z.QPVerified)), ...
        'MaxAbsoluteObjectiveDifferenceCommon',finite_max(z.AbsoluteObjectiveDifference(common)), ...
        'MaxRelativeObjectiveDifferenceCommon',finite_max(z.RelativeObjectiveDifference(common)), ...
        'MaxRelativeStateVectorDifferenceCommon',finite_max(z.RelativeStateVectorDifference(common)), ...
        'MeanTSRelativeStateErrorCommon',finite_mean(z.TSRelativeStateError(common)), ...
        'MeanQPRelativeStateErrorCommon',finite_mean(z.QPRelativeStateError(common)), ...
        'MaxFreshResidualDifference',finite_max([z.TSFreshResidualDifference;z.QPFreshResidualDifference]), ...
        'TSCountChecksPassed',nnz(z.TSCountMatch),'RuntimeRankingReported',false);
    rows=push(rows,row);
end
S=struct2table(rows);
end

function v=nested_value(r,field,key)
v=NaN;if isfield(r,field)&&isfield(r.(field),key),v=r.(field).(key);end
end
function rows=push(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end
function v=finite_max(x)
x=x(isfinite(x));if isempty(x),v=NaN;else,v=max(x);end
end
function v=finite_mean(x)
x=x(isfinite(x));if isempty(x),v=NaN;else,v=mean(x);end
end
function h=sha256_file(p)
fid=fopen(p,'rb');assert(fid>=0,'PMUQP:HashRead','Cannot read %s.',p);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    b=fread(fid,1048576,'*uint8');if ~isempty(b),md.update(typecast(b,'int8'));end
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function restore_session(p,r,t)
path(p);rng(r);maxNumCompThreads(t);
end
