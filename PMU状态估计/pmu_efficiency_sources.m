function sources=pmu_efficiency_sources(systems)
%PMU_EFFICIENCY_SOURCES Select the formal records used by the manuscript.
% Prefer completed current-profile runs. Retained records are selected only
% through the explicit manifest, with parameters, backend and kernels checked.
here=fileparts(mfilename('fullpath'));c=pmu_parameters();
if nargin<1,systems=c.systems;end
if ischar(systems)||isstring(systems),systems=cellstr(systems);end
assert(all(ismember(systems,c.systems)),'PMU:Systems','Unknown test system.');
manifest=jsondecode(fileread(fullfile(here,'package_manifest.json')));
sources=struct('System',{},'Path',{},'RecordedProfile',{},'Retained',{});
for j=1:numel(systems)
 sys=systems{j};prefix=sprintf('final_budget%d_rep%d',c.budget,c.repetitions);
 candidates={fullfile(here,'results',[prefix '_' c.profile '.mat']), ...
  fullfile(here,'results',[prefix '_' strjoin(c.main_systems,'_') '_' c.profile '.mat']), ...
  fullfile(here,'results',[prefix '_' sys '_' c.profile '.mat'])};
 chosen='';
 for k=1:numel(candidates)
  if ~isfile(candidates{k}),continue;end
  d=load(candidates{k},'protocol','rows','completed');
  if ~any(strcmp({d.rows.System},sys)),continue;end
  assert(d.completed,'PMU:IncompleteResults','Formal results are incomplete: %s',candidates{k});
  chosen=candidates{k};break;
 end
 if isempty(chosen)
  assert(isfield(manifest,'ManuscriptEfficiency') ...
   &&isfield(manifest.ManuscriptEfficiency.Sources,sys),'PMU:MissingResults', ...
   'No current formal results or documented retained source for %s.',sys);
  relative=manifest.ManuscriptEfficiency.Sources.(sys);
  assert(startsWith(relative,'results/')&&~contains(relative,'..'),'PMU:SourcePath','Invalid result source.');
  [~,stem]=fileparts(relative);stem=regexprep(stem,'_summary$','');
  chosen=fullfile(here,'results',[stem '.mat']);
  assert(isfile(chosen),'PMU:MissingResults','Formal results not found: %s',chosen);
  d=load(chosen,'protocol','rows','completed');
  assert(d.completed,'PMU:IncompleteResults','Formal results are incomplete: %s',chosen);
 end
 p=d.protocol.Config;
 for f={'lambda','target','feasibility_target','budget','repetitions','initial','names'}
  assert(isequaln(p.(f{1}),c.(f{1})),'PMU:RecordedConfig', ...
   'Recorded %s differs from the current configuration in %s.',f{1},chosen);
 end
 assert(strcmp(backend(p,sys),backend(c,sys)),'PMU:RecordedBackend','Recorded projector differs.');
 for k=1:numel(c.params)
  assert(isequal(fieldnames(c.params{k}),fieldnames(p.params{k})),'PMU:RecordedParameters','Parameter fields differ.');
  for f=fieldnames(c.params{k}).'
   a=c.params{k}.(f{1});b=p.params{k}.(f{1});
   if isa(a,'function_handle'),same=isa(b,'function_handle')&&strcmp(func2str(a),func2str(b));
   else,same=isequaln(a,b);end
   assert(same,'PMU:RecordedParameters','Recorded method %d parameter %s differs.',k,f{1});
  end
 end
 kernels={'TSITTCGA_precision.m','AT2DFP.m','ADCGPM.m','TSIM_DFPM_PRP.m', ...
  'HTTCGP.m','huber_irls_same_model.m','pmu_map.m','pmu_run_cgpm.m'};
 for k=1:numel(kernels)
  ix=find(strcmp({d.protocol.SourceSnapshot.Name},kernels{k}));
  assert(isscalar(ix)&&strcmp(d.protocol.SourceSnapshot(ix).Text,fileread(fullfile(here,kernels{k}))), ...
   'PMU:RecordedSource','Numerical kernel differs from formal results: %s',kernels{k});
 end
 rows=d.rows(strcmp({d.rows.System},sys));
 assert(numel(rows)==45*numel(c.names)*c.repetitions&&numel(unique([rows.Case]))==45, ...
  'PMU:FormalCoverage','Expected 45 cases and all method repetitions for %s.',sys);
 sources(end+1)=struct('System',sys,'Path',chosen,'RecordedProfile',p.profile, ...
  'Retained',~strcmp(p.profile,c.profile)); %#ok<AGROW>
end
end
function value=backend(c,sys)
value=c.backend;
if isfield(c,'backend_by_system')&&isfield(c.backend_by_system,sys),value=c.backend_by_system.(sys);end
end
