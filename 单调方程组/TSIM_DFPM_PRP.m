function [Tcpu,NF,Itr,NormF,info] = TSIM_DFPM_PRP(PR, x0, para)
% TSIM-DFPM-PRP solver. See the package README for source equations and choices.

if nargin < 3 || isempty(para)
    para = struct('Itr_max',2000, 'gamma',0.37, 'sigma',1e-4, 'tau',0.5, ...
                  'rho',1.999, 'quiet',true);
end

timer_id = tic;
x0 = x0(:);
x_return = x0;
failure_reason = 'iteration_limit';
ls_max = getfield_def(para,'ls_max',2000);

k_max     = getfield_def(para,'Itr_max', 2000);
gamma     = getfield_def(para,'gamma',   0.37);
sigma     = getfield_def(para,'sigma',   1e-4);
tau       = getfield_def(para,'tau',     0.5);
rho       = getfield_def(para,'rho',     1.999);
quiet     = getfield_def(para,'quiet',   true);
c = getfield_def(para,'c',0.01);
r = getfield_def(para,'r',0.01);
epsilon_seq=getfield_def(para,'epsilon_seq',@(j)1/j^2);
inertial_fraction = getfield_def(para,'inertial_fraction',0.5);
lambda1=getfield_def(para,'lambda1',0.1); lambda2=getfield_def(para,'lambda2',0.1);
lambda3=getfield_def(para,'lambda3',0.76);
lambda=getfield_def(para,'lambda',0.001); nu=getfield_def(para,'nu',0.8);

epsilon  = 1e-6;
if ~quiet
    fprintf('%s & TSIM-DFPM-PRP (JAMC) & gamma=%.4f & sigma=%.4f & tau=%.4f & rho=%.4f\n', ...
        PR.name, gamma, sigma, tau, rho);
end

Fx0     = PR.f(x0, 1);
NF      = 1;
NormFx0 = norm(Fx0);
x00     = x0;
y1      = x0;      % ϖ_{-1} 占位（原文 Step 0：x_0 = x_{-1} = ϖ_{-1} = ϖ_{-2}）
y2      = x0;      % ϖ_{-2} 占位
L1      = 0;

for k = 1:k_max
    if ~isreal(Fx0) || any(~isfinite(Fx0)) || ~isfinite(norm(Fx0))
        L1=-1; failure_reason='nonfinite_initial_residual'; break;
    end

    if k == 1 && NormFx0 <= epsilon
        L1 = 1; NormF = NormFx0; x_return=x0; break;
    end

    j = k-1;
    if j == 0
        y0 = x0;
    else
        ek = epsilon_seq(j);
        ck = c; rk_inertia = 0;
        if norm(x0-y1)>0, ck = min(c,ek/(j^2*norm(x0-y1))); end
        if j>=2
            rk_inertia = r;
            if norm(x00-y2)>0, rk_inertia = min(r,ek/(j^2*norm(x00-y2))); end
        end
        y0 = x0 + inertial_fraction*ck*(x0-y1) ...
            + inertial_fraction*rk_inertia*(x00-y2);
    end
    Fy0     = PR.f(y0, 1);
    NF      = NF + 1;
    NormFy0 = norm(Fy0);
    if ~isreal(Fy0) || any(~isfinite(Fy0)) || ~isfinite(NormFy0)
        L1=-1; failure_reason='nonfinite_inertial_residual'; break;
    end
    % 停机准则：只在当前迭代点 x_k 判断 ||F(x_k)|| <= epsilon；惯性点 ϖ_k 处不终止。

    if k == 1
        dk = -Fy0;
    else
        w0 = Fy0 - Fy0_old;
        NormFy0_old = norm(Fy0_old);
        betak = (Fy0'*w0)/NormFy0_old^2;
        denominator = max(abs(betak)*Normdk^2, max(norm(w0)^2, NormFy0_old^2));
        vk1 = (Fy0'*dk)/denominator;
        vk2 = (Fy0'*Fy0_old)/denominator;
        vk3 = (Fy0'*w0)/denominator;
        dk  = -Fy0 + lambda1*betak*vk1*dk + lambda2*vk2*Fy0_old - lambda3*vk3*w0;
    end

    Normdk = norm(dk);
    if ~isreal(dk) || any(~isfinite(dk)) || ~isfinite(Normdk)
        L1=-1; failure_reason='nonfinite_direction'; break;
    end
    Normdk2 = Normdk^2;
    Fy0_old = Fy0;

    t=gamma; ls_ok=false;
    for ls_it=0:ls_max
        z0_new=y0+t*dk; Fz0_new=PR.f(z0_new,1); NF=NF+1;
        NormFz0=norm(Fz0_new); cf=max(lambda,min(nu,NormFz0)); gtd=-Fz0_new'*dk;
        ls_ok=accepted_step(Fz0_new,gtd,sigma*t*cf*Normdk2);
        if ls_ok, break; end
        tn=t*tau;
        if tn<=0 || tn==t || isequal(z0_new,y0), break; end
        t=tn;
    end
    if ~ls_ok || isequal(z0_new,y0)
        L1=-1; failure_reason='line_search_failure'; break;
    end

    if NormFz0==0
        z1=y0;                   % F(z_k)=0 -> xi_k=0
    else
        uz=Fz0_new/NormFz0;
        z1=y0-rho*(uz'*(y0-z0_new))*uz;
    end
    x1  = PR.f(z1, 2);              % P_C
    Fx1 = PR.f(x1, 1);
    NF  = NF + 1;
    NormFx1 = norm(Fx1);
    if ~isreal(Fx1) || any(~isfinite(Fx1)) || ~isfinite(NormFx1)
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
if L1==1
    Itr=k; info.status='converged';
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
