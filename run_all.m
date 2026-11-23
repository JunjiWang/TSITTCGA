function run_all(mode)
%RUN_ALL  One entry point for the complete numerical study in this folder.
%
%   run_all            use the frozen Group C instances, run A and B, then
%                      print the saved PMU summaries (no re-solve)
%   run_all('full')    rerun Group C with its default design first, then A, B,
%                      and the PMU solve in final mode
%
% Data sharing:
%   Group C  generates the 200 frozen benchmark instances (results/compare.mat)
%   Group A  reads the same instances
%   Group B  reads the same instances
%   PMU      all five algorithms read the same frozen PMU files in PMU状态估计/data
root=fileparts(mfilename('fullpath'));
saved=pwd; cleanup=onCleanup(@()cd(saved)); %#ok<NASGU>
addpath(root);
if nargin<1, mode='reuse'; end
assert(any(strcmp(mode,{'reuse','full'})),'Use reuse or full.');
t0=tic;
if strcmp(mode,'full')
    fprintf('== Group C: monotone benchmarks (full design) ==\n');
    cd(fullfile(root,'单调方程组')); clear main; main();
    fprintf('Group C done (%.1f s)\n',toc(t0));
else
    assert(isfile(fullfile(root,'单调方程组','results','compare.mat')), ...
        'No Group C instances; run run_all(''full'') first.');
    fprintf('== Group C: reuse existing frozen instances ==\n');
end
fprintf('== Group A: line-search criteria ==\n');
cd(fullfile(root,'A组_线搜索准则')); clear main; main();
fprintf('Group A done (%.1f s)\n',toc(t0));
fprintf('== Group B: step-size strategies ==\n');
cd(fullfile(root,'B组_步长策略')); clear main; main();
fprintf('Group B done (%.1f s)\n',toc(t0));
fprintf('== PMU application ==\n');
pmu=fullfile(root,'PMU状态估计'); cd(pmu);
existing=dir(fullfile(pmu,'results','*_summary.csv'));
if strcmp(mode,'full')
    clear run_pmu; run_pmu('final');
elseif ~isempty(existing)
    fprintf('PMU results present (%d summary files):\n',numel(existing));
    for k=1:min(3,numel(existing))
        T=readtable(fullfile(existing(k).folder,existing(k).name));
        fprintf('  %s\n',existing(k).name);
        disp(T);
    end
else
    clear run_pmu; run_pmu('check');
end
fprintf('PMU done (%.1f s)\n',toc(t0));
fprintf('ALL_DONE\n');
end

