function [Tcpu,NF,Itr,NormF,info] = TSITTCGA_precision(PR, x0, para)
% TSITTCGA —— 论文 Algorithm 1 的实现（PMU 状态估计用的共享核心）。
% 记号对应（论文 -> 代码字段）：
%   sigma -> sigma（线搜索参数）, mu -> mu（(line) 中 max{mu,...} 的下界，not mu_k）,
%   kappa -> kappa（收缩比，也用于 warm start min{1,t_{k-1}/kappa}）,
%   a -> a, rho -> gamma（投影松弛 gamma）, eta -> eta, varsigma -> varsigma,
%   tbar -> tau_bar, c -> c, r -> r, epsilon_seq(k) -> epsilon_k,
%   rho_ls -> varrho（k=0 的初始试探步长）。
% 停机准则：仅在当前迭代点判 ||F(x_k)||<=epsilon；或达到 k_max。
%   惯性点 w_k 处不凭残差停机。近零方向时只核验可行候选 P_C(w_k)；
%   候选未达标且方向非零时，继续原线搜索并正常更新方向历史。

if nargin < 3 || isempty(para)
    para = struct;
end

timer_id = tic;
x0 = x0(:);
x_return = x0;
failure_reason = 'iteration_limit';
ls_max = getfield_def(para,'ls_max',2000);

k_max    = getfield_def(para,'Itr_max', 2000);
sigma    = getfield_def(para,'sigma', 0.8);
rho_ls   = getfield_def(para,'rho_ls', 0.8);
mu       = getfield_def(para,'mu',      .6);
kappa    = getfield_def(para,'kappa',   .5);
a        = getfield_def(para,'a',       .1);
varsigma = getfield_def(para,'varsigma',2);
rho      = getfield_def(para,'rho',     1.99);
eta      = getfield_def(para,'eta',     .001);
tbar     = getfield_def(para,'tbar',    .05);
quiet    = getfield_def(para,'quiet',   true);

c = getfield_def(para,'c',.2);
r = getfield_def(para,'r',.05);
epsilon_seq = getfield_def(para,'epsilon_seq',@(j)1/(j+1)^2);
xbar = getfield_def(para,'xbar',x0); xbar = xbar(:);

if ~quiet
    fprintf('%s & TSITTCGA & sigma=%.3g & rho=%.3g\n',PR.name,sigma,rho);
end

epsilon   = getfield_def(para,'epsilon',1e-6);        % Function-norm tolerance

Fx0     = PR.f(x0, 1);
NF      = 1;
NormFx0 = norm(Fx0);
x00     = x0;
y1      = x0;      % ϖ_{-1} 占位（稿件 k=0,1 只用单步惯性）
y2      = x0;      % ϖ_{-2} 占位
L1      = 0;

t_prev  = rho_ls;
dx_prev = 0;      % ||x_{k-1}-w_{k-2}||，供下一个迭代的 r_k 复用
NormFprev = 0;    % ||F(w_{k-1})||，避免重复计算

debug_dir=logical(getfield_def(para,'debug_dir',false));
if debug_dir, dir_hist=struct('k',{},'j',{},'x',{},'w',{},'F',{},'d',{}); end
for k = 1:k_max
    if k == 1 && (~isreal(Fx0) || ~isfinite(NormFx0))
        L1=-1; failure_reason='nonfinite_initial_residual'; break;
    end

    if k == 1 && NormFx0 <= epsilon
        L1 = 1; NormF = NormFx0; x_return=x0; break;
    end

    j = k-1;
    if j == 0
        y0 = x0 + c*(x0-xbar);
    else
        ek = epsilon_seq(j);
        ck = c; rk_inertia = 0;
        dx = norm(x0-y1);
        if dx>0, ck = min(c,ek/(j^2*dx)); end
        if j>=2
            rk_inertia = r;
            if dx_prev>0, rk_inertia = min(r,ek/(j^2*dx_prev)); end
        end
        y0 = x0 + ck*(x0-y1) + rk_inertia*(x00-y2);
        dx_prev = dx;
    end
    Fy0     = PR.f(y0, 1);
    NF      = NF + 1;
    NormFy0 = norm(Fy0);
    if ~isreal(Fy0) || ~isfinite(NormFy0)
        L1=-1; failure_reason='nonfinite_inertial_residual'; break;
    end

    % Algorithm 1 terminates only at the current iterate x_k.  There is no
    % termination at the inertial point varpi_k: even when ||F(varpi_k)|| is
    % small, the loop continues with the direction, the line search and the
    % projection update below.

    if k == 1
        dk = -Fy0;
    else
        dk_prev     = dk;
        Normd_prev  = Normdk;
        NormFy0_old = NormFprev;
        yk_prev     = Fy0 - Fy0_old;
        sk_prev     = y0  - y0_prev;
        s2          = dot(sk_prev,sk_prev);
        yk_norm_sq  = dot(yk_prev,yk_prev);

        if s2 == 0 || yk_norm_sq == 0
            dk = -Fy0;
        else
            sy    = dot(sk_prev,yk_prev);
            Ahat  = -(3/5)*sy;
            rk    = max(0, eta*NormFy0 - (sy+Ahat)/s2);
            ytil  = yk_prev + (eta+Ahat/s2+rk)*sk_prev;
            nytil = norm(ytil);
            tspec = min(tbar, max(0, 1 - sy/yk_norm_sq));
            mu_k  = max([NormFy0_old^2, varsigma*Normd_prev*nytil, dot(dk_prev,ytil)]);
            if mu_k <= 0
                L1=-1; failure_reason='nonpositive_direction_denominator'; break;
            end
            fdp    = dot(Fy0,dk_prev);
            betak  = dot(Fy0,ytil)/mu_k - (nytil^2*fdp)/mu_k^2;
            thetak = tspec*fdp/mu_k;
            dk     = -Fy0 + betak*dk_prev + thetak*ytil;
        end
    end

    Normdk = norm(dk);
    if ~isreal(dk) || any(~isfinite(dk)) || ~isfinite(Normdk)
        L1=-1; failure_reason='nonfinite_direction'; break;
    end
    if debug_dir
        % 供 verify_algorithm.m 独立复算方向与惯性点（不影响数值结果）。
        dir_hist(end+1)=struct('k',k,'j',k-1,'x',x0,'w',y0,'F',Fy0,'d',dk); %#ok<AGROW>
    end
    Normdk2  = Normdk^2;
    Fy0_old  = Fy0;
    NormFprev = NormFy0;
    y0_prev  = y0;

    if k==1, t=rho_ls; else, t=min(1,t_prev/kappa); end
    z0_new=y0+t*dk; Fz0_new=PR.f(z0_new,1); NF=NF+1;
    NormFz0=norm(Fz0_new);
    cf=max(norm(Fz0_new-Fy0)/(a+NormFz0+NormFy0),mu);
    gtd=-dot(Fz0_new,dk);
    ls_ok=accepted_step(Fz0_new,NormFz0,gtd,sigma*t*cf*Normdk2);
    if ls_ok
        if t<1 && NormFz0>0
            t_hat=min(1,t/kappa);
            z0_hat=y0+t_hat*dk; Fz0_hat=PR.f(z0_hat,1); NF=NF+1;
            NormFzhat=norm(Fz0_hat);
            cfh=max(norm(Fz0_hat-Fy0)/(a+NormFzhat+NormFy0),mu);
            gtdh=-Fz0_hat'*dk;
            if accepted_step(Fz0_hat,NormFzhat,gtdh,sigma*t_hat*cfh*Normdk2)
                t=t_hat; z0_new=z0_hat; Fz0_new=Fz0_hat; NormFz0=NormFzhat;
            end
        end
    else
        for ls_it=1:ls_max
            tn=kappa*t;
            if tn<=0 || tn==t || isequal(z0_new,y0), break; end
            t=tn; z0_new=y0+t*dk; Fz0_new=PR.f(z0_new,1); NF=NF+1;
            NormFz0=norm(Fz0_new);
            cf=max(norm(Fz0_new-Fy0)/(a+NormFz0+NormFy0),mu);
            gtd=-dot(Fz0_new,dk);
            ls_ok=accepted_step(Fz0_new,NormFz0,gtd,sigma*t*cf*Normdk2);
            if ls_ok, break; end
        end
    end
    if ~ls_ok || isequal(z0_new,y0)
        L1=-1; failure_reason='line_search_failure'; break;
    end
    t_prev=t;

    if NormFz0==0
        z1=y0;                   % F(z_k)=0 -> xi_k=0, x_{k+1}=P_C(varpi_k)
    else
        proj_coef=rho*dot(Fz0_new,y0-z0_new)/(NormFz0^2);
        z1=y0-proj_coef*Fz0_new;
    end
    x1  = PR.f(z1, 2);              % P_C
    Fx1 = PR.f(x1, 1);
    NF  = NF + 1;
    NormFx1 = norm(Fx1);
    if ~isreal(Fx1) || ~isfinite(NormFx1)
        L1=-1; failure_reason='nonfinite_projected_residual'; break;
    end
    if NormFx1 <= epsilon
        L1 = 1; NormF = NormFx1; x_return=x1; break;
    end

    x00 = x0;  x0 = x1;  NormFx0 = NormFx1; Fx0=Fx1;
    y2  = y1;  y1 = y0;          % ϖ_{k-2} ← ϖ_{k-1} ← ϖ_k
end

Tcpu=toc(timer_id);
actual_NF=NF; actual_time=Tcpu;
if L1==1
    % 停机只看残差；可行性只记录，不参与成功/失败的判定。
    feasibility=norm(x_return-PR.f(x_return,2));
    if ~isreal(x_return) || any(~isfinite(x_return)) || ~isfinite(NormF) ...
            || NormF>epsilon
        L1=-1; failure_reason='invalid_termination';
    end
else
    feasibility=NaN; x_return=x0;
end
info=struct('success',L1==1,'status',failure_reason,'actual_NF',actual_NF, ...
    'actual_time',actual_time,'loop_iterations',k,'feasibility',feasibility,'x',x_return);
if debug_dir, info.dir_history=dir_hist; end
if L1==1
    Itr=k; info.status='converged';
else
    NF=NaN; Itr=NaN; Tcpu=NaN; NormF=NaN;
end
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function ok=accepted_step(Fz,fz,lhs,rhs)
% A feasible exact root terminates the solve before a projection update.
% fz = ||Fz|| is passed in: the caller already computed it, and a finite norm
% implies that every component of Fz is finite.
ok=isreal(Fz) && isfinite(fz) && isreal(lhs) && isfinite(lhs) ...
    && isreal(rhs) && isfinite(rhs) && lhs>=rhs;
end
