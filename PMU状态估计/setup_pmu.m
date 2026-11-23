function setup_pmu(silent)
%SETUP_PMU  One-time path setup for the PMU state-estimation package.
%   setup_pmu          clean stale entries, add this folder, and save the path
%   setup_pmu(true)    same, with less console output
%
% Why: the ABC benchmark folders of the old workspace (TSITTCGA_A_B_C数值实验)
% contain files with the same names as the ones used here (pmu_parameters,
% pmu_prepare_backend, TSITTCGA, main, ...).  While those folders stay on the
% MATLAB search path they shadow this package and "main" either is not found
% or runs the wrong code.  This script removes them, puts the package first
% on the path, and saves the result so the next MATLAB session keeps it.
if nargin<1, silent=false; end
here=fileparts(mfilename('fullpath'));
p=path;
parts=strsplit(p,pathsep);
keep=true(size(parts));
stale={'TSITTCGA_A_B_C数值实验','PMU状态估计_替换前备份'};
for k=1:numel(parts)
    for s=1:numel(stale)
        if contains(parts{k},stale{s}) && ~strcmpi(parts{k},here)
            keep(k)=false;break
        end
    end
end
removed=parts(~keep);
if any(~keep)
    rmpath(strjoin(removed,pathsep));
end
addpath(here,'-begin');
addpath(fullfile(here,'tools'),'-begin');
try
    savepath;
    saved=true;
catch err
    saved=false;
    fprintf('Path could not be saved (%s). Add the folder manually each session.\n',err.message);
end
if ~silent
    fprintf('Package folder : %s\n',here);
    fprintf('Removed stale path entries: %d\n',numel(removed));
    for k=1:numel(removed), fprintf('  - %s\n',removed{k}); end
    if saved, fprintf('Search path saved for future MATLAB sessions.\n'); end
    fprintf('\nReady. Type pmu_main for the menu, or main for a full rerun.\n\n');
end
end
