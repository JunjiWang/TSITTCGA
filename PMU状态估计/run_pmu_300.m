function T=run_pmu_300(mode,repetitions)
%RUN_PMU_300  IEEE300-only run with exactly the same frozen configuration.
%   run_pmu_300            final mode, three timed repetitions
%   run_pmu_300('check')   the three check cases only
% The parameter set, lambda, tolerance, budget and initial point are taken
% from pmu_parameters.m, so this run is directly comparable with the
% IEEE14/IEEE118 main program.
here=fileparts(mfilename('fullpath'));
if nargin<1||isempty(mode), mode='final'; end
if nargin<2, repetitions=[]; end
c=pmu_parameters();
T=run_pmu(mode,c.large_system,c.budget,repetitions);
end
