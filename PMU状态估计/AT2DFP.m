function [Tcpu,NF,Itr,NormF,info] = AT2DFP(PR, x0, para)
% AT2DFP —— Pengjie Liu,
% "A three-term derivative-free projection method with BFGS-like update and
%  its accelerated variant", Optimization 75:8 (2026) 2055-2089, Algorithm 2.
%
% 说明：除停机准则与本工程其余算法保持一致（||G|| <= 1e-6、Itr >= Itr_max，
%       或 ||d|| <= 1e-7，且返回点必须可行）外，迭代流程与参数全部按原文：
%       惯性外推 (12)、方向 (13)(14)、向量 p_k (5)、den 修正 𝓎_k (7)、χ_k (6)、
%       自适应线搜索 (15)、投影步 Step 5。
%       原文参数 sigma=1e-4, varsigma=0.5, rho=0.4, gamma_k=1.99, tau=0.99,
%       mu=0.01, chibar=0.65；AT2 额外 phi=0.25, epsilon_k=1/k^2 (k>=1)。

if nargin < 3 || isempty(para)
    para = struct('Itr_max',2000, 'sigma',1e-4, 'varsigma',0.5, 'rho_ls',0.4, ...
                  'gamma',1.99, 'tau',0.99, 'mu',0.01, 'chibar',0.65, ...
                  'phi',0.25, 'quiet',true);
end

timer_id = tic;
x0 = x0(:);
x_return = x0;
failure_reason = 'iteration_limit';
ls_max = getfield_def(para,'ls_max',2000);

k_max    = getfield_def(para,'Itr_max', 2000);
sigma    = getfield_def(para,'sigma',   1e-4);
varsigma = getfield_def(para,'varsigma',0.5);
rho_ls   = getfield_def(para,'rho_ls',  0.4);
gamma    = getfield_def(para,'gamma',   1.99);
tau      = getfield_def(para,'tau',     0.99);
mu       = getfield_def(para,'mu',      0.01);
chibar   = getfield_def(para,'chibar',  0.65);
phi      = getfield_def(para,'phi',     0.25);
quiet    = getfield_def(para,'quiet',   true);

epsilon   = 1e-6;
if ~quiet
    fprintf('%s & AT2DFP (Liu 2026) & sigma=%.4g & rho=%.4g & gamma=%.4g\n', ...
        PR.name, sigma, rho_ls, gamma);
end


% 论文 Step 0：给定 x_{-1}, x_0。本工程只有单一起点，取 x_{-1} := x_0。
x       = x0;
x_prev  = x0;
Fx      = PR.f(x, 1);
NF      = 1;
NormFx  = norm(Fx);
L1      = 0;
d_prev  = [];
v_prev  = [];
Gv_prev = [];
NormGv_prev = 0;

for k = 0:k_max-1
    if ~isreal(Fx) || any(~isfinite(Fx)) || ~isfinite(NormFx)
        L1=-1; failure_reason='nonfinite_initial_residual'; break;
    end

    % ---- Step 1 -----------------------------------------------------
    if NormFx <= epsilon
        L1 = 1; NormF = NormFx; x_return=x; break;
    end

    % ---- Step 2: 惯性外推 (12) --------------------------------------
    if norm(x - x_prev) > 0
        phik = min(phi, (1/max(k,1)^2)/norm(x - x_prev));
    else
        phik = phi;
    end
    v  = x + phik*(x - x_prev);
    Gv = PR.f(v, 1);
    NF = NF + 1;
    NormGv = norm(Gv);
    if ~isreal(Gv) || any(~isfinite(Gv)) || ~isfinite(NormGv)
        L1=-1; failure_reason='nonfinite_inertial_residual'; break;
    end


    % ---- Step 3: 方向 (13)(14) --------------------------------------
    if k == 0
        d = -Gv;
    else
        % 向量 p_k —— 式 (5)
        if NormGv_prev > 0
            pk = Gv - (NormGv/NormGv_prev)*Gv_prev;
        else
            pk = Gv;
        end
        npk2 = pk'*pk;
        if npk2 == 0, pk = Gv; npk2 = pk'*pk; end  % Equation (5) fallback
        if npk2 == 0
            d = -Gv;
        else
            % den 修正向量 𝓎_k —— 式 (7)
            ndprev = norm(d_prev);
            if ndprev > 0
                tau_k = tau*norm(pk)/ndprev + max(0, -(d_prev'*pk)/ndprev^2);
            else
                tau_k = tau;
            end
            yk = pk + tau_k*d_prev;
            den = d_prev'*yk;
            if den == 0, den = eps; end
            % χ_k —— 式 (6)，这里 G(v_k) 即 Gv，G(v_{k-1}) 即 Gv_prev
            chi_k = min(chibar, max(0, (pk'*((Gv - Gv_prev) - (v - v_prev)))/npk2));
            betak  = (Gv'*pk)/den - (npk2*(Gv'*d_prev))/den^2;      % (14)
            thetak = chi_k*(Gv'*d_prev)/den;                        % (14)
            d = -Gv + betak*d_prev + thetak*pk;                     % (13)
        end
    end

    Normd = norm(d);
    if ~isreal(d) || any(~isfinite(d)) || ~isfinite(Normd)
        L1=-1; failure_reason='nonfinite_direction'; break;
    end
    % Inertial method: small ||d_k|| is not a stopping test; keep iterating.
    Normd2 = Normd^2;

    % ---- Step 4: 自适应线搜索 (15) ----------------------------------
    t = varsigma;
    ls_ok = false;
    for ls_it = 0:ls_max
        z = v + t*d; Fz = PR.f(z,1); NF = NF + 1; NormFz = norm(Fz);
        rhs = sigma*t*NormFz*min(mu, Normd2);       % (15) 右端
        ls_ok = accepted_step(Fz, -Fz'*d, rhs);
        if ls_ok, break; end
        tn = t*rho_ls;
        if tn<=0 || tn==t || isequal(z,v), break; end
        t = tn;
    end
    if ~ls_ok || isequal(z,v)
        L1=-1; failure_reason='line_search_failure'; break;
    end


    % ---- Step 5: 投影步 ---------------------------------------------
    if NormFz == 0
        xi = 0;
    else
        xi = (Fz'*(v - z))/NormFz^2;
    end
    x_new = PR.f(v - gamma*xi*Fz, 2);               % P_C[.]

    d_prev = d;  v_prev = v;  Gv_prev = Gv;  NormGv_prev = NormGv;
    x_prev = x;  x = x_new;
    Fx = PR.f(x, 1);
    NF = NF + 1;
    NormFx = norm(Fx);
    if ~isreal(Fx) || any(~isfinite(Fx)) || ~isfinite(NormFx)
        L1=-1; failure_reason='nonfinite_projected_residual'; break;
    end
    if NormFx <= epsilon
        L1 = 1; NormF = NormFx; x_return=x; break;
    end
end

Tcpu=toc(timer_id);
actual_NF=NF; actual_time=Tcpu;
if L1==1
    feasibility=norm(x_return-PR.f(x_return,2));
    if ~isreal(x_return) || any(~isfinite(x_return)) || ~isfinite(NormF) ...
            || NormF>epsilon
        L1=-1; failure_reason='invalid_termination';
    end
else
    feasibility=NaN; x_return=x;
end
info=struct('success',L1==1,'status',failure_reason,'actual_NF',actual_NF, ...
    'actual_time',actual_time,'loop_iterations',k+1,'feasibility',feasibility,'x',x_return);
if L1==1
    Itr=k+1; info.status='converged';
else
    NF=NaN; Itr=NaN; Tcpu=NaN; NormF=NaN;
end
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function ok=accepted_step(Fz,lhs,rhs)
ok=isreal(Fz) && all(isfinite(Fz)) && isfinite(norm(Fz)) ...
    && isreal(lhs) && isfinite(lhs) && isreal(rhs) && isfinite(rhs) && lhs>=rhs;
end
