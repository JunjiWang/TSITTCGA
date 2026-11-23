function main(action)
%MAIN  Rerun the PMU state-estimation experiments from scratch.
%
%   main            archive the old results and rerun everything:
%                   IEEE14 + IEEE118, IEEE300, WLS accuracy, QP check,
%                   figures and the exported tables.  Several hours.
%   main('quick')   quick acceptance run only (IEEE14 all cases,
%                   IEEE118/300 cases 1,16,31).  About six minutes.
%
% Old result files are moved to results/archive/full_<stamp>/ first.  If a
% file cannot be moved, another MATLAB experiment is probably still running
% and the run stops with a clear message instead of mixing results.
here=fileparts(mfilename('fullpath'));
addpath(here,'-begin');
if nargin<1||isempty(action), action='full'; end
switch lower(action)
    case {'full','all'}
        run_full(here);
    case 'quick'
        archive_results(here,{'check_*'});
        run_pmu('check');
    otherwise
        error('main:Action','Unknown action "%s". Use main or main(''quick'').',action);
end
end

function run_full(here)
t0=tic;
fprintf('=== PMU experiments: full rerun ===\n');
fprintf('Started: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
archive_results(here,{'*'});

fprintf('[1/6] IEEE14 + IEEE118, six methods ...\n');
run_pmu('final');

fprintf('\n[2/6] IEEE300, same parameter set ...\n');
run_pmu_300();

fprintf('\n[3/6] WLS paired accuracy ...\n');
run_pmu_wls_comparison();

fprintf('\n[4/6] QP reference check ...\n');
run_pmu_qp_validation();

fprintf('\n[5/6] figures ...\n');
pmu_figures('main');
try, pmu_figures('300'); catch err, fprintf('  (300 profiles skipped: %s)\n',err.message); end
pmu_state_figures('ieee14');

fprintf('\n[6/6] exported tables ...\n');
addpath(fullfile(here,'tools'));
try
    write_pmu_supplement_tables();
catch err
    fprintf('  (tables skipped: %s)\n',err.message);
end

fprintf('\n=== finished in %.1f minutes ===\n',toc(t0)/60);
fprintf('Results : %s\n',fullfile(here,'results'));
fprintf('Figures : %s\n',fullfile(here,'figures'));
end

function archive_results(here,patterns)
out=fullfile(here,'results');
dst=fullfile(out,'archive',['full_' datestr(now,'yyyymmdd_HHMMSS')]);
files=[];
for k=1:numel(patterns)
    f=dir(fullfile(out,patterns{k}));
    f=f(~[f.isdir]);
    files=[files; f]; %#ok<AGROW>
end
if isempty(files)
    fprintf('No old result files to archive.\n\n');
    return
end
if ~isfolder(dst), mkdir(dst); end
failed={};
for k=1:numel(files)
    src=fullfile(out,files(k).name);
    [ok,msg]=movefile(src,dst,'f');
    if ~ok
        failed{end+1}=sprintf('%s  (%s)',files(k).name,msg); %#ok<AGROW>
    end
end
if ~isempty(failed)
    error('main:Archive', ...
        ['Could not archive %d of %d result file(s):\n  %s\n\n' ...
         'Another MATLAB experiment is probably still writing them.\n' ...
         'Stop that run (Ctrl-C in its window), delete or move the files, ' ...
         'and call main again.'], numel(failed), numel(files), strjoin(failed,newline));
end
fprintf('Archived %d old result file(s) to\n  %s\n\n',numel(files),dst);
end
