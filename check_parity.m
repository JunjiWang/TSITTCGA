function check_parity()
% Verify that the controlled solver copy matches the native one bitwise.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'单调方程组'),'-begin');
addpath(fullfile(root,'线搜索对比'),'-begin');
s=load(fullfile(root,'单调方程组','results','compare.mat'),'inst','X0','wjj');
p=s.wjj;
cfg=struct('criterion','Proposed','sigma',p.sigma,'initial',p.rho_ls, ...
    'shrink',p.kappa,'warm',true,'increase',true,'lambda',.001,'nu',.8, ...
    'liu_mu',.01,'j_a',2,'j_l',.5,'label','Proposed');
nrep=size(s.X0,2); nc=numel(s.inst)*nrep;
ok=true;
for i=[1,nc]
    q=ceil(i/nrep); r=i-nrep*(q-1);
    short=p; short.Itr_max=8;
    [~,~,~,~,a]=TSITTCGA(s.inst(q),s.X0{q,r},short);
    [~,~,~,~,b]=TSITTCGA_LS(s.inst(q),s.X0{q,r},short,cfg);
    same=isequaln(a.x,b.x) && a.actual_NF==b.actual_NF ...
        && a.loop_iterations==b.loop_iterations && a.success==b.success;
    fprintf('case %d (inst %d, start %d): same=%d xdiff=%.3e NF %d/%d itr %d/%d\n', ...
        i,q,r,same,norm(a.x-b.x),a.actual_NF,b.actual_NF,a.loop_iterations,b.loop_iterations);
    ok=ok&&same;
end
fprintf('PARITY_%d\n',ok);
end
