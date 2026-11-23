function report=verify_pmu_results(profile)
% Independently recompute saved endpoints using fresh sparse least squares.
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');
c=pmu_parameters();
current_sources=nargin<1;
if current_sources,profile=c.profile;end
data=cell(1,3);
for j=1:3
    z=load(fullfile(here,'data',['input_' c.systems{j} '.mat']),'instances');data{j}=z.instances;
end
if current_sources
    selected=pmu_efficiency_sources();
    paths=unique({selected.Path},'stable');files=struct('name',{});
    for j=1:numel(paths),[~,stem,ext]=fileparts(paths{j});files(j).name=[stem ext];end
else
    files=[dir(fullfile(out,['check_budget*_' profile '.mat'])); ...
        dir(fullfile(out,['final_budget*_' profile '.mat']))];
end
assert(~isempty(files),'No completed comparison records found.');
entries=struct([]);pending={};total=0;cgchecks=0;successes=0;retained_failures=0;maxdiff=0;maxobj=0;max_success=0;
for f=1:numel(files)
    d=load(fullfile(out,files(f).name),'rows','endpoints','completed','protocol');
    if ~d.completed
        pending{end+1}=files(f).name; %#ok<AGROW>
        fprintf('Pending, excluded from endpoint audit: %s\n',files(f).name);
        continue
    end
    assert(numel(d.rows)==numel(d.endpoints),'Endpoint count mismatch.');
    recorded=d.protocol.Config;
    if current_sources
        selected_systems={selected(strcmp({selected.Path},fullfile(out,files(f).name))).System};
    else
        assert(strcmp(recorded.profile,profile),'Selected file has a different profile.');
        selected_systems=c.systems;
    end
    verified_file_runs=0;
    for k=1:numel(d.rows)
        if ~ismember(d.rows(k).System,selected_systems),continue;end
        verified_file_runs=verified_file_runs+1;
        row=d.rows(k);r=d.endpoints{k};sy=find(strcmp(c.systems,row.System));
        inst=data{sy}{row.Case};m=inst.m;s=r.s(:);e=s(1:m)-s(m+1:end);
        xi=sparse(inst.G)\(inst.y-e);fit=inst.y-e-inst.G*xi;
        F=[min(s(1:m),recorded.lambda-fit);min(s(m+1:end),recorded.lambda+fit)];
        fn=norm(F);feas=norm(min(s,0));obj=.5*dot(fit,fit)+recorded.lambda*norm(e,1);
        delta=abs(fn-row.Residual);objdelta=abs(obj-row.Objective);
        assert(all(isfinite([s;xi;F]))&&delta<=1e-9&&objdelta<=1e-9,'Independent endpoint mismatch.');
        assert(norm(e-r.e)==0,'Saved bad-data estimate differs from returned point.');
        assert(row.Success==r.success,'Saved success flags differ.');
        if row.Success
            assert(fn<=recorded.target&&feas<=recorded.feasibility_target,'Successful endpoint fails original criteria.');
            max_success=max(max_success,fn);
            successes=successes+1;
        else,retained_failures=retained_failures+1;end
        if ~strcmp(row.Method,'Huber-IRLS')
            assert(r.count_ok&&r.report_ok&&r.NFCounted==r.info.actual_NF,'NF or reported residual mismatch.');
            cgchecks=cgchecks+1;
        end
        total=total+1;maxdiff=max(maxdiff,delta);maxobj=max(maxobj,objdelta);
    end
    entry=struct('File',files(f).name,'Runs',verified_file_runs,'Completed',d.completed,'RecordedProfile',recorded.profile);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
end
assert(total>0,'PMU:NoCompletedResults','No completed runs for the selected profile.');
a=load(fullfile(out,'data_audit.mat'),'report');assert(a.report.Summary.Pass);
report=struct('Created',char(datetime('now')),'Pass',true,'InputCases',a.report.Summary.PassedCases, ...
    'ParameterProfile',profile,'PendingFiles',{pending},'AllSelectedFilesCompleted',isempty(pending), ...
    'SavedRunsVerified',total,'SuccessfulRuns',successes,'RetainedFailureRuns',retained_failures, ...
    'NativeCGPMCountChecks',cgchecks,'MaxFreshResidualDifference',maxdiff, ...
    'MaxFreshObjectiveDifference',maxobj,'MaxFreshSuccessfulResidual',max_success,'Files',entries, ...
    'Scope','Selected formal records verified against frozen input data; retained source profiles are listed per file.');
save(fullfile(out,'endpoint_audit.mat'),'report');
fid=fopen(fullfile(out,'endpoint_audit.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));disp(report);
end
