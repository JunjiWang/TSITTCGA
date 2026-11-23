function main(options)
% Group A. Default reuses verified existing results; resume=false reruns.
% main(struct('figs_only',true)) redraws without solving.
% main(struct('criterion_cost',false)) skips the criterion cost study.
if nargin==0, options=struct; end
cost_study=~isfield(options,'criterion_cost') || options.criterion_cost;
if isfield(options,'criterion_cost'), options=rmfield(options,'criterion_cost'); end
assert(~isfield(options,'mode'),'Use the B-group entry for initial-step experiments.');
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
    assert(isfile(fullfile(here,'results','compare.mat')), ...
        'GroupA:MissingResults','Run Group A before requesting figures only.');
else
    assert(isfile(fullfile(root,'单调方程组','results','compare.mat')), ...
        'GroupA:MissingGroupC','Run 单调方程组/main.m first to create the shared test instances.');
end
options.phases={'A'};
main_line_search(options);
if cost_study && ~(isfield(options,'figs_only') && options.figs_only)
    criterion_cost(struct('source_root',root,'outdir',fullfile(here,'results')));
end
end
