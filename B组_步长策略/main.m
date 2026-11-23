function main(options)
% Group B compares fixed-initial reduction with warm start plus one increase test.
% main(struct('figs_only',true)) redraws without solving.
% main(struct('mode','initial_steps')) runs the whole-solve initial-step ablation.
% main(struct('mode','step_probe')) compares matched first-search branches.
if nargin==0, options=struct; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
shared=fullfile(root,'线搜索对比');
savedpath=path; cleanup=onCleanup(@()path(savedpath)); 
addpath(shared,'-begin');
assert(strcmpi(which('main_line_search'),fullfile(shared,'main_line_search.m')), ...
    'Unexpected experiment driver on the MATLAB path.');
assert(~isfield(options,'phases') && ~isfield(options,'output_dir'), ...
    'This entry point fixes its group and output folder.');
options.output_dir=here;
if ~isfield(options,'resume'), options.resume=true; end
if isfield(options,'figs_only') && options.figs_only
    result_name='compare.mat';
    if isfield(options,'mode') && any(strcmp(options.mode,{'initial_steps','step_probe'}))
        result_name='initial_steps.mat';
    end
    assert(isfile(fullfile(here,'results',result_name)), ...
        'GroupB:MissingResults','Run the requested Group B experiment before requesting figures only.');
elseif ~isfield(options,'mode') || ~strcmp(options.mode,'step_probe')
    assert(isfile(fullfile(root,'单调方程组','results','compare.mat')), ...
        'GroupB:MissingGroupC','Run 单调方程组/main.m first to create the shared test instances.');
end
if isfield(options,'mode') && any(strcmp(options.mode,{'initial_steps','step_probe'}))
    main_line_search(options); return
end
options.phases={'B'};
main_line_search(options);
end
