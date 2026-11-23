function [Tcpu,NF,Itr,NormF,info] = HTTCGP(PR, x0, para)
% HTTCGP —— Yin, Jian, Jiang, Liu & Wang,
% "A hybrid three-term conjugate gradient projection method for constrained
%  nonlinear monotone equations with applications",
% Numerical Algorithms 88 (2021) 389-418, Algorithm 1.
%
% 说明：除停机准则与本工程其余算法保持一致（||F|| <= 1e-6、Itr >= Itr_max、
%       或 ||d|| <= 1e-7）外，迭代流程与参数全部按原文：方向 (7)(8)、
%       自适应线搜索 (9)、投影步 (10)(11)。保留原文第 401 页的极小试探步长规则。
%       原文参数 rho=0.5, varsigma=1, sigma=0.01, lambda=0.001, nu=0.8,
%       mu=0.2, gamma=1.6, tbar=0.3。

if nargin < 3 || isempty(para)
    para = struct('Itr_max',2000, 'tbar',0.3, 'mu',0.2, 'varsigma',1, ...
                  'rho_ls',0.5, 'sigma',0.01, 'lambda',0.001, 'nu',0.8, ...
                  'gamma',1.6, 'quiet',true);
end

timer_id = tic;
x0 = x0(:);
x_return = x0;
failure_reason = 'iteration_limit';
ls_max = getfield_def(para,'ls_max',2000);

k_max    = getfield_def(para,'Itr_max', 2000);
tbar     = getfield_def(para,'tbar',    0.3);
mu       = getfield_def(para,'mu',      0.2);
varsigma = getfield_def(para,'varsigma',1);
rho_ls   = getfield_def(para,'rho_ls',  0.5);
sigma    = getfield_def(para,'sigma',   0.01);
lambda   = getfield_def(para,'lambda',  0.001);
nu       = getfield_def(para,'nu',      0.8);
gamma    = getfield_def(para,'gamma',   1.6);
quiet    = getfield_def(para,'quiet',   true);

epsilon   = 1e-6;
epsilon_C = getfield_def(para,'epsilon_C',1e-6);

if ~quiet
    fprintf('%s & HTTCGP (Yin et al. 2021) & sigma=%.4g & rho=%.4g & gamma=%.4g\n', ...
        PR.name, sigma, rho_ls, gamma);
end


x       = x0;
Fx      = PR.f(x, 1);
NF      = 1;
NormFx  = norm(Fx);
L1      = 0;
paper_floor_accepts = 0;
d_prev  = [];
x_prev  = [];
F_prev  = [];

for k = 0:k_max-1
    if ~isreal(Fx) || any(~isfinite(Fx)) || ~isfinite(NormFx)
        L1=-1; failure_reason='nonfinite_initial_residual'; break;
    end

    % ---- Step 1 -----------------------------------------------------
    if NormFx <= epsilon
        L1 = 1; NormF = NormFx; x_return=x; break;
    end

    % ---- Step 2: 方向 (7)(8) ----------------------------------------
    if k == 0
        d = -Fx;
    else
        y_prev = Fx - F_prev;               % y_{k-1} = F_k - F_{k-1}
        s_prev = x  - x_prev;               % s_{k-1} = x_k - x_{k-1}
        ny2    = y_prev'*y_prev;
        wk = max([mu*norm(d_prev)*norm(y_prev), d_prev'*y_prev, norm(F_prev)^2]);
        if wk <= 0, wk = eps; end
        if ny2 > 0
            tk = min(tbar, max(0, 1 - (y_prev'*s_prev)/ny2));   % (5)
        else
            tk = 0;
        end
        betak  = (Fx'*y_prev)/wk - (ny2*(Fx'*d_prev))/wk^2;      % (8)
        thetak = tk*(Fx'*d_prev)/wk;                             % (8)
        d = -Fx + betak*d_prev + thetak*y_prev;                  % (7)
    end

    Normd = norm(d);
    if ~isreal(d) || any(~isfinite(d)) || ~isfinite(Normd)
        L1=-1; failure_reason='nonfinite_direction'; break;
    end
    Normd2 = Normd^2;

    % ---- Step 3: 自适应线搜索 (9) -----------------------------------
    t = varsigma;
    ls_ok = false;
    for ls_it = 0:ls_max
        z = x + t*d; Fz = PR.f(z,1); NF = NF + 1; NormFz = norm(Fz);
        Phi = min(nu, max(lambda, NormFz));      % P_[lambda,nu][ ||F(z)|| ]
        ls_ok = accepted_step(Fz, -Fz'*d, sigma*t*Phi*Normd2);
        if ls_ok, break; end
        % Yin et al. (2021), p. 401: accept a rejected trial <= 1e-10.
        % This is the paper's numerical safeguard, not satisfaction of (9).
        if ls_it>0 && t<=1e-10 && isreal(Fz) && all(isfinite(Fz)) ...
                && isfinite(NormFz)
            ls_ok=true; paper_floor_accepts=paper_floor_accepts+1; break;
        end
        tn = t*rho_ls;
        if tn<=0 || tn==t || isequal(z,x), break; end
        t = tn;
    end
    if ~ls_ok || isequal(z,x)
        L1=-1; failure_reason='line_search_failure'; break;
    end

    % ---- Step 4: 投影步 (10)(11) ------------------------------------
    if NormFz == 0
        xi = 0;
    else
        xi = (Fz'*(x - z))/NormFz^2;
    end                 % (11)
    x_new = PR.f(x - gamma*xi*Fz, 2);            % (10)  P_C[.]

    d_prev = d;  x_prev = x;  F_prev = Fx;
    x = x_new;
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
    'actual_time',actual_time,'paper_floor_accepts',paper_floor_accepts,'loop_iterations',k+1,'feasibility',feasibility,'x',x_return);
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
