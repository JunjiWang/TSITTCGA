function result=criterion_cost(options)
% Same-data timing of the proposed and Yin (2021) line-search criteria.
% All trials of the 40 standard problem/dimension pairs are recorded
% without altering native solvers. Low-dimensional problems and signal
% recovery are excluded.
if nargin==0, options=struct; end
code_dir=fileparts(mfilename('fullpath')); here=code_dir;
default_root='C:\Users\Administrator\Desktop\TSITTCGA_A_B_C数值实验';
if isfolder(fullfile(fileparts(code_dir),'单调方程组')), default_root=fileparts(code_dir); end
o=struct('source_root',default_root,'outdir','','repeats',12,'target_seconds',.025, ...
    'threads',1,'case_indices',[],'save_results',true,'plots',true,'refresh_only',false);
for field=fieldnames(options)'
    assert(isfield(o,field{1}),'Unknown option: %s',field{1});
    o.(field{1})=options.(field{1});
end
assert(o.repeats>=4 && mod(o.repeats,2)==0 && o.target_seconds>0);
if ~isempty(o.outdir), here=o.outdir; if ~isfolder(here), mkdir(here); end, end
if o.refresh_only
    data=load(fullfile(here,'results.mat'),'result'); result=data.result;
    write_outputs(result,here,o.plots); return
end
savedpath=path; savedrng=rng; savedthreads=maxNumCompThreads;
cleanup=onCleanup(@()restore_session(savedpath,savedrng,savedthreads)); %#ok<NASGU>
maxNumCompThreads(o.threads); addpath(code_dir,'-begin');
root=o.source_root; groups={'单调方程组'};
native_sources=cellfun(@(g)fileread(fullfile(root,g,'TSITTCGA.m')),groups,'UniformOutput',false);
probe_sources=cellfun(@make_probe_source,native_sources,'UniformOutput',false);
[cases,input_meta]=make_cases(root,groups);
existing_A=read_existing_A(root);
if ~isempty(o.case_indices), cases=cases(o.case_indices); end
source_files={}; source_snapshot={};
for g=1:numel(groups)
    for file={'main.m','TSITTCGA.m','problems.m','problem.m'}
        f=fullfile(root,groups{g},file{1});
        if isfile(f), source_files{end+1}=f; source_snapshot{end+1}=fileread(f); end %#ok<AGROW>
    end
end
names={'Proposed','Yin2021'};
scenarios={'Cached norms: adaptive term','Cached values: criterion','Function evaluation + criterion'};
rows=cell(numel(cases),24); raw=cell(numel(cases),1); probe_group=0;
wall=tic;
fprintf('CURRENT CRITERION TIMING: %d cases, %d paired blocks, %d thread(s).\n',numel(cases),o.repeats,o.threads);
for i=1:numel(cases)
    c=cases(i); folder=fullfile(root,groups{c.group});
    if c.group~=probe_group
        clear probe
        probe=load_cost_probe(code_dir,probe_sources{c.group}); probe_group=c.group;
    end
    addpath(folder,'-begin'); clear TSITTCGA
    native=str2func('TSITTCGA');
    assert(strcmpi(functions(native).file,fullfile(folder,'TSITTCGA.m')),'Wrong native solver on path.');
    [~,nf,it,fn,d]=native(c.PR,c.x0,c.p); expected=signature(nf,it,fn,d);
    assert(d.success,'The prespecified case failed; do not silently drop it.');
    bench=struct('mode','record','values',[],'samples',true,'capture_all',true);
    [~,nfp,itp,fnp,dp]=probe(c.PR,c.x0,c.p,bench);
    assert(isequaln(expected,signature(nfp,itp,fnp,dp)),'Recording changed the native numerical result.');
    packets=dp.cost_samples; nt=numel(packets);
    assert(nt==dp.diff_count && nt==numel(dp.diff_values) && nt>0);
    assert(isequal(cellfun(@(s)s.index,packets),1:nt));
    ratios=zeros(1,nt);
    for q=1:nt
        s=packets{q}; u=c.PR.f(s.z,1);
        assert(isequaln(u,s.u) && isequal(norm(u-s.v),s.difference_norm));
        assert(all(isfinite(u)) && all(isfinite(s.v)));
        s.nz=norm(u); packets{q}=s;
        ratios(q)=s.difference_norm/(c.p.a+s.nz+s.nv);
    end
    % The same actual trials and cached base norm are used in both criteria.
    reference=zeros(2,3);
    for scenario=1:3
        for kind=1:2
            [reference(kind,scenario),calls]=cost_batch(c.PR,packets,c.p,kind,scenario,1);
            assert(calls==(scenario==3)*nt);
            cost_batch(c.PR,packets,c.p,kind,scenario,1); % untimed warm-up
        end
    end
    seconds=nan(o.repeats,2,3); loops=zeros(1,3);
    for scenario=1:3
        pilot=zeros(1,2);
        for kind=1:2
            timer=tic; cost_batch(c.PR,packets,c.p,kind,scenario,1);
            pilot(kind)=toc(timer);
        end
        % Calibrate to the slower criterion. Using the faster scalar clip
        % would repeat the high-dimensional difference norm unnecessarily.
        % Both criteria still execute exactly the same number of passes.
        loops(scenario)=min(20000,max(1,ceil(o.target_seconds/max(pilot))));
        for rep=1:o.repeats
            order=1:2; if mod(i+rep+scenario,2), order=2:-1:1; end
            for kind=order
                timer=tic;
                [value,calls]=cost_batch(c.PR,packets,c.p,kind,scenario,loops(scenario));
                seconds(rep,kind,scenario)=toc(timer)/loops(scenario);
                target=reference(kind,scenario)*loops(scenario);
                assert(abs(value-target)<=1e-9*max(1,abs(target)),'Timed arithmetic checksum mismatch.');
                assert(calls==(scenario==3)*nt*loops(scenario),'Timed function count mismatch.');
            end
        end
    end
    % Independent uninstrumented solver timing supplies a cost reference.
    timer=tic; native_batch(native,c,1); pilot=toc(timer);
    native_loops=min(500,max(1,ceil(o.target_seconds/max(pilot,eps))));
    native_seconds=zeros(o.repeats,1);
    for rep=1:o.repeats
        timer=tic; actual=native_batch(native,c,native_loops);
        native_seconds(rep)=toc(timer)/native_loops;
        assert(isequaln(expected,actual),'A native timed run changed its output.');
    end
    med=squeeze(median(seconds,1));
    extra=seconds(:,1,1)-seconds(:,2,1);
    criterion_extra=seconds(:,1,2)-seconds(:,2,2);
    trial_pct=100*(seconds(:,1,3)-seconds(:,2,3))./seconds(:,2,3);
    pct=100*extra./native_seconds;
    iteration=cellfun(@(s)s.iteration,packets); steps=cellfun(@(s)s.t,packets);
    initial=[true,diff(iteration)~=0]; increase=[false,diff(iteration)==0 & diff(steps)>0];
    rows(i,:)={c.label,c.group_name,numel(c.x0),c.start,it,nf,nt,sum(initial), ...
        nt-sum(initial)-sum(increase),sum(increase),median(native_seconds), ...
        med(1,1)/nt*1e6,med(2,1)/nt*1e6,median(extra)/nt*1e6, ...
        med(1,2)/nt*1e6,med(2,2)/nt*1e6,median(criterion_extra)/nt*1e6, ...
        med(1,3)/nt*1e6,med(2,3)/nt*1e6,median(trial_pct), ...
        median(pct),local_quantile(pct,.25),local_quantile(pct,.75),true};
    raw{i}=struct('case',c.label,'group',c.group_name,'n',numel(c.x0), ...
        'timings',seconds,'native_seconds',native_seconds,'batch_loops',loops, ...
        'native_batch_loops',native_loops,'reference_checksums',reference, ...
        'signature',{expected},'parameters',parameter_snapshot(c.p),'x0',c.x0, ...
        'trial_steps',steps,'trial_iterations',iteration,'difference_values',dp.diff_values, ...
        'cached_base_norms',cellfun(@(s)s.nv,packets), ...
        'cached_trial_norms',cellfun(@(s)s.nz,packets), ...
        'source_sha256',sha256(native_sources{c.group}),'probe_sha256',sha256(probe_sources{c.group}), ...
        'all_trials_recorded',true,'native_exact_match',true);
    fprintf('%2d/%d %-10s n=%6d LS=%4d factor extra=%8.3fus estimated/native=%+.2f%% [%.2f,%.2f] elapsed=%.1fs\n', ...
        i,numel(cases),c.label,numel(c.x0),nt,median(extra)/nt*1e6,median(pct), ...
        local_quantile(pct,.25),local_quantile(pct,.75),toc(wall));
    if o.save_results
        progress=cell2table(rows(1:i,:),'VariableNames',cost_columns);
        writetable(progress,fullfile(here,'case_costs.csv'),'Encoding','UTF-8');
    end
    clear packets dp
end
for j=1:numel(source_files)
    assert(strcmp(source_snapshot{j},fileread(source_files{j})),'A formal source changed during timing.');
end
if existing_A.available
    assert(strcmp(existing_A.sha256,hash_file(existing_A.file)),'Existing Group A results changed during timing.');
end
result=struct('cases',cell2table(rows,'VariableNames',cost_columns),'raw',{raw}, ...
    'options',o,'criteria',{names},'scenarios',{scenarios},'input_meta',input_meta, ...
    'source_files',{source_files},'source_snapshot',{source_snapshot}, ...
    'probe_sources',{probe_sources},'driver_source',fileread([mfilename('fullpath'),'.m']), ...
    'existing_A',existing_A, ...
    'date',char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss')), ...
    'matlab_version',version,'cpu',getenv('PROCESSOR_IDENTIFIER'),'elapsed_seconds',toc(wall), ...
    'all_numerical_checks_passed',true);
if o.save_results
    save(fullfile(here,'results.mat'),'result','-v7'); write_outputs(result,here,o.plots);
end
fprintf('DONE: all %d original solver/recorded solver comparisons matched exactly.\n',numel(cases));
disp(dimension_summary(result.cases));
end

function names=cost_columns
names={'Case','Group','n','Start','Itr','NF','LSNF','InitialTrials','BacktrackingTrials', ...
    'IncreaseTests','NativeSeconds','ProposedFactorUs','YinFactorUs','FactorExtraUs', ...
    'ProposedCriterionUs','YinCriterionUs','CriterionExtraUs', ...
    'ProposedFunctionCriterionUs','YinFunctionCriterionUs','FunctionCriterionExtraPct', ...
    'EstimatedExtraPctOfNative','ExtraPctQ25','ExtraPctQ75','NativeExactMatch'};
end

function [value,calls]=cost_batch(PR,packets,p,kind,scenario,loops)
value=0; calls=0;
for j=1:loops
    for q=1:numel(packets)
        s=packets{q}; u=s.u; nz=s.nz;
        if scenario==3, u=PR.f(s.z,1); calls=calls+1; end
        if scenario>=2, nz=norm(u); end
        if kind==1
            factor=max(p.mu,norm(u-s.v)/(p.a+nz+s.nv));
        else
            factor=max(.001,min(.8,nz));
        end
        if scenario==1
            value=value+factor;
        else
            lhs=-dot(u,s.d); rhs=p.sigma*s.t*factor*s.nd2;
            value=value+factor+double(lhs>=rhs);
        end
    end
end
end

function T=dimension_summary(S)
keys=unique(S(:,{'Group','n'}),'rows','stable'); rows=cell(height(keys),9);
for j=1:height(keys)
    take=strcmp(S.Group,keys.Group{j}) & S.n==keys.n(j); a=S(take,:);
    rows(j,:)={keys.Group{j},keys.n(j),height(a),sum(a.LSNF), ...
        median(a.ProposedFactorUs),median(a.YinFactorUs),median(a.FactorExtraUs), ...
        median(a.CriterionExtraUs),median(a.EstimatedExtraPctOfNative)};
end
T=cell2table(rows,'VariableNames',{'Group','n','Cases','Trials','ProposedFactorUs', ...
    'YinFactorUs','FactorExtraUs','CriterionExtraUs','EstimatedExtraPctOfNative'});
end

function write_outputs(r,here,plots)
S=r.cases; D=dimension_summary(S);
writetable(D,fullfile(here,'dimension_summary.csv'),'Encoding','UTF-8');
fid=fopen(fullfile(here,'实验说明与结果.md'),'w','n','UTF-8'); assert(fid>0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# 线搜索准则计算开销实验\n\n');
fprintf(fid,'状态：VERIFIED（数值轨迹与计数核对）；计时结果为本机实测估计。\n\n日期：%s，Asia/Shanghai。MATLAB：%s，线程数：%d。\n\n',r.date,r.matlab_version,r.options.threads);
fprintf(fid,'## 实验范围\n\n复用当前正式代码、参数及保存的初始点。40个标准问题与维数组合各取原保存的第一个起点，共40例，不包含低维问题和信号恢复。若使用case_indices运行子集，以实际表格为准。没有重新调参。所有原生求解器、问题文件、主程序与TeX保持原样。\n\n');
fprintf(fid,'每例记录TSITTCGA的全部实际线搜索试探，涵盖初始试探、缩步试探和单次步长增大检验。记录版本逐项核对原生版本的返回点、Itr、NF、最终函数范数、状态与实际求值次数，所有计时中的重复求解也检查相同输出。\n\n');
fprintf(fid,'## 三种计时范围\n\n');
fprintf(fid,'1. 自适应项：复用已计算的F(z)、F(w)、两点范数，比较max{mu, norm(F(z)-F(w))/(a+norm(F(z))+norm(F(w)))}与Yin的max{lambda,min{nu,norm(F(z))}}。Yin取lambda=0.001、nu=0.8。差向量及其范数在每次计时中实时计算，没有预存替代。\n');
fprintf(fid,'2. 准则判断：复用函数值，双方均重新计算试探点范数、内积与右端表达式；基点范数与方向范数均缓存。\n');
fprintf(fid,'3. 函数求值加准则判断：双方均在相同的记录试探点重新计算F(z)，再执行第2种计算。每块实际求值次数逐项核对一致。\n\n');
fprintf(fid,'每种场景采用%d个配对计时块，两准则先后顺序均衡。先预热，再按同一批次数重复整段试探序列；根据较慢项校准至约%.3f秒，较快项的计时块可能更短。校验计算结果的累加值，防止计时未执行目标运算。原生求解器另外重复计时%d次。保留全部重复、负时间差及四分位范围，不截取有利结果。\n\n',r.options.repeats,r.options.target_seconds,r.options.repeats);
fprintf(fid,'## 按维数汇总\n\n时间单位为微秒，每例先在重复计时中取中位数，再在同组同维数算例间取中位数。新增时间先逐个配对块相减，再取中位数，因此不必等于另外两栏中位数的差。\n\n');
fprintf(fid,'| 组别 | 映射维数 | 算例数 | 本文自适应项(us) | Yin自适应项(us) | 单次新增(us) | 累计新增/原生求解估计 |\n|---|---:|---:|---:|---:|---:|---:|\n');
for j=1:height(D)
    fprintf(fid,'| %s | %d | %d | %.4f | %.4f | %+.4f | %+.3f%% |\n',D.Group{j},D.n(j),D.Cases(j),D.ProposedFactorUs(j),D.YinFactorUs(j),D.FactorExtraUs(j),D.EstimatedExtraPctOfNative(j));
end
fprintf(fid,'\n## 按实验组汇总\n\n');
fprintf(fid,'| 组别 | 算例数 | 每次新增(us)中位数 | 累计新增/原生求解中位数 | 函数求值加判断相对Yin增幅中位数 |\n|---|---:|---:|---:|---:|\n');
for group=unique(S.Group,'stable')'
    take=strcmp(S.Group,group{1}); a=S(take,:);
    fprintf(fid,'| %s | %d | %+.4f | %+.3f%% | %+.3f%% |\n',group{1},height(a),median(a.FactorExtraUs),median(a.EstimatedExtraPctOfNative),median(a.FunctionCriterionExtraPct));
end
fprintf(fid,'\n累计新增时间以相同全部试探序列中两种自适应项计算时间的配对差估计，除以另外测得的原生完整求解时间。它衡量相同数据上的额外算术成本，不能解释为将正式线搜索替换为Yin后算法总Time的变化。记录试探的集中重复访问会改变缓存行为，因此该比例是成本估计，不能视为无误差的分项时间。\n\n');
fprintf(fid,'## 完整求解中的净效果\n\n现有A组比较保持相同算法、方向、初始试探步长及缩减因子，关闭热启动与增大检验。以下直接读取当前A组结果，不重新选择参数或删除算例。\n\n');
if r.existing_A.available
    a=r.existing_A;
    fprintf(fid,'| 设置 | 成功数 | Itr均值 | NF均值 | LSNF均值 | Time均值(s) |\n|---|---:|---:|---:|---:|---:|\n');
    for j=1:2
        fprintf(fid,'| %s | %d/%d | %.4f | %.4f | %.4f | %.6f |\n',r.criteria{j},a.success(j),a.ncase,a.Itr(j),a.NF(j),a.LSNF(j),a.Time(j));
    end
    fprintf(fid,'\nA组原结果文件SHA-256：`%s`，此摘要已冻结保存于本实验results.mat中。\n',a.sha256);
end
fprintf(fid,'\nA组的Time来自其原来整组实验，本次微观计时来自另一次运行，不能跨批次相减。应分别报告每次判断的额外计算成本及A组完整求解表现。LSNF只统计映射求值，不增加向量加减、范数或标量运算的计数。\n\n');
fprintf(fid,'## 当前参数与复现\n\n');
group_names={'Monotone'};
for g=1:numel(group_names)
    ix=find(cellfun(@(q)strcmp(q.group,group_names{g}),r.raw),1);
    if isempty(ix), continue; end
    p=r.raw{ix}.parameters;
    fprintf(fid,'- %s：sigma=%.6g, mu=%.6g, a=%.6g, kappa=%.6g, 初始步长=%.6g, varsigma=%.6g, gamma=%.6g, eta=%.6g, bar_tau=%.6g, c=%.6g, r=%.6g。源码SHA-256：`%s`。\n',r.raw{ix}.group,p.sigma,p.mu,p.a,p.kappa,p.rho_ls,p.varsigma,p.rho,p.eta,p.tbar,p.c,p.r,r.raw{ix}.source_sha256);
end
fprintf(fid,'\n在此文件夹运行main。results.mat保存全部配对计时、参数、起点、源码快照和核对结果。case_costs.csv为逐例结果，dimension_summary.csv为按维数汇总。TSITTCGA_cost_probe.m是从当前原生源码自动生成的记录助手，不是另一种对比算法。\n');
if plots, plot_costs(D,here); end
end

function plot_costs(D,here)
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17.6 6.4]);
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile; hold(ax,'on'); take=strcmp(D.Group,'Monotone'); a=D(take,:);
plot(ax,a.n,a.ProposedFactorUs,'-o','Color',[0 .40 .70],'MarkerSize',4,'LineWidth',1.4,'DisplayName','Proposed');
plot(ax,a.n,a.YinFactorUs,'--s','Color',[.10 .60 .45],'MarkerSize',4,'LineWidth',1.4,'DisplayName','Yin (2021)');
xlabel(ax,'Dimension'); ylabel(ax,'Time per adaptive term (\mus)');
title(ax,'(a) Standard equations'); legend(ax,'Location','northwest','Box','off','FontSize',8);
grid(ax,'on'); box(ax,'on');
ax=nexttile; hold(ax,'on'); colors=[0 .40 .70;.85 .33 .10;.10 .60 .45];
groups={'Monotone'}; labels={'Standard equations'};
for g=1:numel(groups)
    a=D(strcmp(D.Group,groups{g}),:);
    scatter(ax,a.n,a.EstimatedExtraPctOfNative,28,colors(g,:),'filled','DisplayName',labels{g});
end
set(ax,'XScale','log'); yline(ax,0,':','HandleVisibility','off');
xlabel(ax,'Mapping dimension'); ylabel(ax,'Estimated extra time / total time (%)');
title(ax,'(b) Recorded-trial cost estimate'); legend(ax,'Location','best','Box','off','FontSize',7);
grid(ax,'on'); box(ax,'on');
set(findall(f,'Type','axes'),'FontName','Times New Roman','FontSize',8);
exportgraphics(f,fullfile(here,'cost_summary.png'),'Resolution',300);
exportgraphics(f,fullfile(here,'cost_summary.pdf'),'ContentType','vector');
end

function a=read_existing_A(root)
file=fullfile(root,'A组_线搜索准则','results','compare.mat');
a=struct('available',false,'file',file);
if ~isfile(file), return; end
d=load(file,'experiments'); assert(isfield(d,'experiments') && isfield(d.experiments,'A'));
s=d.experiments.A; a.available=true; a.sha256=hash_file(file); a.ncase=size(s.NF,1);
a.labels=s.labels(1:2); a.success=sum(s.Success(:,1:2),1);
for key={'Itr','NF','LSNF','Time'}
    a.(key{1})=mean(s.(key{1})(:,1:2),1);
end
end

function h=hash_file(file)
fid=fopen(file,'rb'); assert(fid>0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1048576,'*uint8'); md.update(bytes);
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function sig=signature(nf,it,f,d)
sig={nf,it,f,d.success,d.status,d.actual_NF,d.loop_iterations,d.x};
end

function sig=native_batch(solver,c,loops)
for j=1:loops
    [~,nf,it,f,d]=solver(c.PR,c.x0,c.p);
end
sig=signature(nf,it,f,d);
end

function [cases,meta]=make_cases(root,groups)
cases=struct('label',{},'group',{},'group_name',{},'start',{},'PR',{},'x0',{},'p',{});
savedpath=path; clean=onCleanup(@()path(savedpath)); %#ok<NASGU>
meta=struct;
for g=1:numel(groups)
    folder=fullfile(root,groups{g}); addpath(folder,'-begin');
    text=fileread(fullfile(folder,'main.m'));
    key='wjj';
    expr=regexp(text,[key,'=struct\([\s\S]*?\);'],'match','once');
    assert(~isempty(expr),'Cannot find current parameter block.');
    expr=regexprep(expr,['^',key,'='],'p=');
    o=struct('Itr_max',2000); %#ok<NASGU> parameters use o.Itr_max
    eval(expr);
    assert(p.Itr_max==2000 && ~isfield(p,'line_search'));
    meta.(['group',num2str(g)])=struct('source',folder,'parameter_expression',expr, ...
        'solver_sha256',sha256(fileread(fullfile(folder,'TSITTCGA.m'))));
    if g==1
        data=load(fullfile(folder,'results','compare.mat'),'inst','X0');
        clear problems; pf=str2func('problems');
        assert(strcmpi(functions(pf).file,fullfile(folder,'problems.m')));
        for j=1:numel(data.inst)
            id=data.inst(j).id; x=data.X0{j,1};
            PR=struct('name',sprintf('P%d',id),'n',numel(x),'f',@(v,m)pf(v,m,id));
            cases(end+1)=make_case(PR,x,p,g,'Monotone',1); %#ok<AGROW>
        end
        assert(numel(data.inst)==40,'Expected all 40 problem-dimension pairs.');
    end
end
end
function c=make_case(PR,x,p,g,name,start)
c=struct('label',PR.name,'group',g,'group_name',name,'start',start, ...
    'PR',PR,'x0',x,'p',p);
end
function source=make_probe_source(native)
source=strrep(native,'= TSITTCGA(PR, x0, para)', ...
    '= TSITTCGA_cost_probe(PR, x0, para, bench)');
assert(~strcmp(source,native),'Unexpected native function signature.');
init=sprintf(['cost_mode=find(strcmp(bench.mode,{''record'',''compute'',''replay''}));\n', ...
    'assert(isscalar(cost_mode)); cost_count=0; cost_values=bench.values;\n', ...
    'cost_samples=cell(1,3); cost_capture=bench.samples;\n', ...
    'cost_capture_all=isfield(bench,''capture_all'') && bench.capture_all;\n']);
source=strrep(source,'timer_id = tic;',[init,'timer_id = tic;']);
assert(numel(strfind(source,'norm(Fz0_new-Fy0)'))==2);
assert(numel(strfind(source,'norm(Fz0_hat-Fy0)'))==1);
source=strrep(source,'norm(Fz0_new-Fy0)', ...
    'cost_difference(Fz0_new,Fy0,z0_new,dk,Normdk2,t)');
source=strrep(source,'norm(Fz0_hat-Fy0)', ...
    'cost_difference(Fz0_hat,Fy0,z0_hat,dk,Normdk2,t_hat)');
footer=sprintf([ ...
    'info.diff_count=cost_count;\n', ...
    'if cost_mode==1, info.diff_values=cost_values; end\n', ...
    'if cost_capture\n', ...
    '    cost_samples=cost_samples(~cellfun(@isempty,cost_samples));\n', ...
    '    [~,cost_unique]=unique(cellfun(@(q)q.index,cost_samples),''stable'');\n', ...
    '    info.cost_samples=cost_samples(cost_unique);\n', ...
    'end\n', ...
    '    function v=cost_difference(u,w,z,d,nd2,t)\n', ...
    '        cost_count=cost_count+1;\n', ...
    '        switch cost_mode\n', ...
    '            case 1\n', ...
    '                v=norm(u-w); cost_values(cost_count)=v;\n', ...
    '            case 2\n', ...
    '                v=norm(u-w);\n', ...
    '            case 3\n', ...
    '                v=cost_values(cost_count);\n', ...
    '        end\n', ...
    '        if cost_capture\n', ...
    '            q=struct(''u'',u,''v'',w,''z'',z,''d'',d,''nd2'',nd2, ...\n', ...
    '                ''t'',t,''nv'',norm(w),''index'',cost_count, ...\n', ...
    '                ''iteration'',k,''difference_norm'',v);\n', ...
    '            if cost_capture_all\n', ...
    '                cost_samples{cost_count}=q;\n', ...
    '            else\n', ...
    '                if cost_count==1, cost_samples{1}=q; end\n', ...
    '                if 2^floor(log2(cost_count))==cost_count, cost_samples{2}=q; end\n', ...
    '                cost_samples{3}=q;\n', ...
    '            end\n', ...
    '        end\n', ...
    '    end\n']);
marker=regexp(source,'(?m)^end\r?\n\r?\nfunction v = getfield_def','start','once');
assert(~isempty(marker));
source=[source(1:marker-1),footer,source(marker:end)];
end
function probe=load_cost_probe(here,source)
% Regenerate from the matching group source; only one temporary probe is kept.
probe_file=fullfile(here,'TSITTCGA_cost_probe.m');
if ~isfile(probe_file) || ~strcmp(fileread(probe_file),source)
    fid=fopen(probe_file,'w','n','UTF-8'); assert(fid>0);
    fprintf(fid,'%s',source); fclose(fid);
end
clear TSITTCGA_cost_probe
rehash;
probe=str2func('TSITTCGA_cost_probe');
assert(strcmpi(functions(probe).file,probe_file),'Wrong probe on MATLAB path.');
end

function q=local_quantile(x,p)
x=sort(x(:)); u=1+(numel(x)-1)*p; lo=floor(u); hi=ceil(u);
q=x(lo)'.*(hi-u)+x(hi)'.*(u-lo); same=lo==hi; q(same)=x(lo(same));
end
function s=parameter_snapshot(s)
for field=fieldnames(s)'
    if isa(s.(field{1}),'function_handle')
        s.(field{1})=func2str(s.(field{1}));
    end
end
end

function h=sha256(text)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(unicode2native(text,'UTF-8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function restore_session(p,r,t)
path(p); rng(r); maxNumCompThreads(t);
end




