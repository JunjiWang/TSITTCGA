function pmu_main(action)
%PMU_MAIN Single entry point for the PMU robust state-estimation package.
%
%   pmu_main                 show the command menu and the file map
%   pmu_main('final')        main comparison: IEEE14 + IEEE118, six methods, 3 reps
%   pmu_main('300')          IEEE300 alone, same frozen parameter set
%   pmu_main('check')        quick acceptance run (IEEE14 all, 118/300 cases 1,16,31)
%   pmu_main('wls')          paired WLS / TSITTCGA / Huber-IRLS accuracy study
%   pmu_main('qp')           independent QP accuracy reference on 10 fixed cases
%   pmu_main('verify')       independent recomputation of the current profile
%   pmu_main('audit')        audit the frozen input data
%
% Every experiment refuses to overwrite an existing result file. Rename or move
% the old file first. Verification never changes the experiment results.
here=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(action)
    show_menu(here);
    return
end
oldpath=path; cleanup=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(here,'-begin');
switch lower(action)
    case 'final',  run_pmu('final');
    case 'check',  run_pmu('check');
    case '300',    run_pmu_300();
    case 'wls',    run_pmu_wls_comparison();
    case 'qp',     run_pmu_qp_validation();
    case 'audit',  audit_pmu_data();
    case 'tables',   addpath(fullfile(here,'tools')); write_pmu_supplement_tables();
    case 'build300'
        addpath(fullfile(here,'tools'));
        try
            build_ieee300_input();
        catch err
            fprintf('\n%s\n',err.message);
            fprintf(['The frozen IEEE300 input in data/ is still valid; ' ...
                     'this tool is only needed to regenerate it.\n']);
        end
    case 'verify', verify_all(here);
    otherwise, error('PMU:Entry','Unknown action "%s". Run pmu_main with no argument for the menu.',action);
end
end

function show_menu(here)
out=fullfile(here,'results');
fprintf('PMU 鲁棒状态估计实验包\n');
fprintf('  数据: %s\n',fullfile(here,'data'));
fprintf('  结果: %s\n',out);
fprintf('\n命令\n');
fprintf('  pmu_main(''final'')    主比较 (IEEE14 + IEEE118, 6 方法, 3 次重复)\n');
fprintf('  pmu_main(''check'')    快速验收 (14 全量, 118/300 取 1,16,31 例)\n');
fprintf('  pmu_main(''300'')      IEEE300 单独运行（同一组参数）\n');
fprintf('  pmu_main(''wls'')      配对精度: WLS / TSITTCGA / Huber-IRLS\n');
fprintf('  pmu_main(''qp'')       独立 QP 精度参照 (固定 10 例)\n');
fprintf('  pmu_main(''verify'')   独立重算并核对当前配置的结果\n');
fprintf('  pmu_main(''audit'')    审计冻结量测数据\n');
fprintf('  pmu_main(''tables'')    导出 WLS/QP 表格到 docs/（不参与实验）\n');
fprintf('  pmu_main(''build300'')  从 MATPOWER 重建 IEEE300 输入（不参与实验）\n');
fprintf('\n命令对应脚本\n');
fprintf('  final/check -> run_pmu.m\n');
fprintf('  wls         -> run_pmu_wls_comparison.m\n');
fprintf('  qp          -> run_pmu_qp_validation.m + pmu_qp_reference.m\n');
fprintf('  verify      -> verify_pmu_results.m / verify_pmu_supplement_results.m\n');
fprintf('  audit       -> audit_pmu_data.m + pmu_build.m\n');
fprintf('\n求解器与模型\n');
fprintf('  TSITTCGA_precision.m, ADCGPM.m, TSIM_DFPM_PRP.m, HTTCGP.m, AT2DFP.m\n');
fprintf('  pmu_run_cgpm.m, pmu_map.m, huber_irls_same_model.m\n');
fprintf('  pmu_prepare_backend.m, pmu_parameters.m, summarize_pmu.m\n');
fprintf('\n文档生成 (不参与实验)\n');
fprintf('  tools/write_pmu_supplement_tables.m, tools/build_ieee300_input.m\n');
end

function verify_all(here)
out=fullfile(here,'results');c=pmu_parameters();n=0;
formal=dir(fullfile(out,['final_budget*_' c.profile '.mat']));
if ~isempty(formal)
    if ~isfile(fullfile(out,'data_audit.mat'))
        data_check=audit_pmu_data();
        assert(data_check.Summary.Pass,'PMU:DataAudit','Frozen input audit failed.');
    end
    r=verify_pmu_results();n=n+1;
    fprintf('[%d] 论文所用正式结果核对: %d 次运行\n',n,r.SavedRunsVerified);
else
    fprintf('主结果未核对：缺少配置 %s 的正式结果。\n',c.profile);
end
if isfile(fullfile(out,['wls_paired_' c.profile '.mat'])) && isfile(fullfile(out,['qp_validation_' c.profile '.mat']))
    verify_pmu_supplement_results();n=n+1;
    fprintf('[%d] 补充实验核对完成 (%s)\n',n,c.profile);
else
    fprintf('补充实验未核对：配置 %s 的 WLS 或 QP 结果缺失。\n',c.profile);
end

if n==0, error('PMU:Entry','No saved results to verify.'); end
fprintf('已完成上述 %d 组核对。\n',n);
end

