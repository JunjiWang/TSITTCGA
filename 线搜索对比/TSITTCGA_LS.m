function [Tcpu,NF,Itr,NormF,info] = TSITTCGA_LS(PR, x0, para, ls)
% Controlled LS experiments. The direction, inertia and projection are
% copied from the current TSITTCGA.m; only the line search is generalized.

if nargin < 3 || isempty(para)
    para = struct;
end

assert(nargin==4 && isstruct(ls),'Explicit LS configuration required.');
assert(any(strcmp(ls.criterion,{'Proposed','Yin2021','Liu','Jiang'})));
assert(ls.sigma>0 && ls.sigma<1 && ls.initial>0 && ls.initial<=1);
assert(ls.shrink>0 && ls.shrink<1 && isscalar(ls.warm) && isscalar(ls.increase));
assert(ls.lambda>0 && ls.nu>=ls.lambda && ls.liu_mu>0 && ls.j_a>1 && ls.j_l>0);
timer_id = tic;
x0 = x0(:);
x_return = x0;
failure_reason = 'iteration_limit';
ls_max = getfield_def(para,'ls_max',2000);
assert(isscalar(ls_max) && isfinite(ls_max) && ls_max>=0 && ls_max==fix(ls_max), 'Invalid ls_max.');

k_max    = getfield_def(para,'Itr_max', 2000);
sigma=ls.sigma;
rho_ls=ls.initial;
kappa=ls.shrink;
mu       = getfield_def(para,'mu',      .6);
a        = getfield_def(para,'a',       .01);
para.mu=mu; para.a=a; % Pass resolved defaults to line_test.
varsigma = getfield_def(para,'varsigma',2);
rho      = getfield_def(para,'rho',     1.98);
eta      = getfield_def(para,'eta',     .01);
tbar     = getfield_def(para,'tbar',    .05);
quiet    = getfield_def(para,'quiet',   true);

c = getfield_def(para,'c',0.05);
r = getfield_def(para,'r',.001);
epsilon_seq = getfield_def(para,'epsilon_seq',@(j)1/(j+1)^2);
xbar = getfield_def(para,'xbar',x0); xbar = xbar(:);
assert(c>0 && c<1 && r>0 && r<1 && numel(xbar)==numel(x0), 'Invalid inertial parameters.');

if ~quiet
    label=ls.label;
    fprintf('%s & %s & sigma=%.3g & rho=%.3g\n',PR.name,label,sigma,rho);
end

epsilon   = 1e-6;                                   % Function-norm tolerance
epsilon_C = getfield_def(para,'epsilon_C', 1e-6);

assert(k_max>=1 && k_max==fix(k_max) && isfinite(k_max), 'Invalid iteration limit.');
assert(isreal(x0) && all(isfinite(x0)) && epsilon_C>0 && isfinite(epsilon_C), 'Invalid initial point/tolerance.');
assert(norm(x0-PR.f(x0,2))<=epsilon_C,'Initial point must be feasible.');
assert(isreal(xbar) && all(isfinite(xbar)) && norm(xbar-PR.f(xbar,2))<=epsilon_C,'xbar must be feasible.');
assert(mu>0 && mu<=1 && a>0 && isfinite(a),'Invalid proposed-criterion parameters.');
assert(rho>0 && rho<2,'Projection relaxation must be in (0,2).');
assert(varsigma>0 && isfinite(varsigma) && eta>0 && isfinite(eta) && tbar>=0 && tbar<1,'Invalid direction parameters.');
Fx0     = PR.f(x0, 1);
NF      = 1;
NormFx0 = norm(Fx0);
x00     = x0;
y1      = x0;      % ϖ_{-1} 占位（稿件 k=0,1 只用单步惯性）
y2      = x0;      % ϖ_{-2} 占位
L1      = 0;
t_prev  = rho_ls;
ls_nf=0; backtracks=0; increase_tests=0; increase_accepts=0;
initial_accepts=0; inertia_nf=0; projection_nf=0;
accepted_steps=nan(1,k_max);

record_history=isfield(PR,'signal_n') && isfield(PR,'true_x');
if record_history, hm=nan(1,k_max); ht=hm; end
for k = 1:k_max
    if ~isreal(Fx0) || any(~isfinite(Fx0)) || ~isfinite(norm(Fx0))
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
        assert(isscalar(ek) && isfinite(ek) && ek>0 && ek<1, 'epsilon_j must be in (0,1).');
        ck = c; rk_inertia = 0;
        if norm(x0-y1)>0, ck = min(c,ek/(j^2*norm(x0-y1))); end
        if j>=2
            rk_inertia = r;
            if norm(x00-y2)>0, rk_inertia = min(r,ek/(j^2*norm(x00-y2))); end
        end
        y0 = x0 + ck*(x0-y1) + rk_inertia*(x00-y2);
    end
    Fy0     = PR.f(y0, 1);
    NF      = NF + 1;
    inertia_nf=inertia_nf+1;
    NormFy0 = norm(Fy0);
    if ~isreal(Fy0) || any(~isfinite(Fy0)) || ~isfinite(NormFy0)
        L1=-1; failure_reason='nonfinite_inertial_residual'; break;
    end

    % Algorithm 1 terminates only at the current iterate x_k.  There is no
    % termination at the inertial point varpi_k.

    if k == 1
        dk = -Fy0;
    else
        dk_prev     = dk;
        Normd_prev  = Normdk;
        NormFy0_old = norm(Fy0_old);
        yk_prev     = Fy0 - Fy0_old;
        sk_prev     = y0  - y0_prev;
        s2          = dot(sk_prev,sk_prev);
        yk_norm_sq  = dot(yk_prev,yk_prev);

        if s2 == 0 || yk_norm_sq == 0
            dk = -Fy0;
        else
            Ahat = -(3/5)*(dot(sk_prev,yk_prev));
            rk   = max(0, eta*NormFy0 - (dot(sk_prev,yk_prev)+Ahat)/s2);
            ytil = yk_prev + (eta+Ahat/s2+rk)*sk_prev;
            nytil = norm(ytil);
            tspec = min(tbar, max(0, 1 - (dot(yk_prev,sk_prev))/yk_norm_sq));
            mu_k = max([NormFy0_old^2, varsigma*Normd_prev*nytil, dot(dk_prev,ytil)]);
            if mu_k <= 0
                L1=-1; failure_reason='nonpositive_direction_denominator'; break;
            end
            betak  = (dot(Fy0,ytil))/mu_k - (nytil^2*(dot(Fy0,dk_prev)))/mu_k^2;
            thetak = tspec*(dot(Fy0,dk_prev))/mu_k;
            dk     = -Fy0 + betak*dk_prev + thetak*ytil;
        end
    end

    Normdk = norm(dk);
    if ~isreal(dk) || any(~isfinite(dk)) || ~isfinite(Normdk)
        L1=-1; failure_reason='nonfinite_direction'; break;
    end
    if Normdk <= 1e-8
        % Degenerate direction: Algorithm 1 then takes xi_k = 0, that is
        % x_{k+1} = P_C(varpi_k).  Execute that step instead of exiting.
        x1 = PR.f(y0,2); Fx1 = PR.f(x1,1); NF = NF + 1; projection_nf = projection_nf + 1; NormFx1 = norm(Fx1);
        if NormFx1 <= epsilon
            L1 = 1; NormF = NormFx1; x_return = x1; break;
        end
        x00 = x0; x0 = x1; NormFx0 = NormFx1; Fx0 = Fx1;
        y2 = y1; y1 = y0;
        continue;
    end
    Normdk2  = Normdk^2;
    Fy0_old  = Fy0;
    y0_prev  = y0;

    if ls.warm && k>1, t=min(1,t_prev/kappa); else, t=rho_ls; end
    z0_new=y0+t*dk; Fz0_new=PR.f(z0_new,1); NF=NF+1; ls_nf=ls_nf+1;
    ls_ok=line_test(PR,z0_new,Fz0_new,Fy0,NormFy0,dk,Normdk2,t,para,ls,epsilon_C);
    if ls_ok
        initial_accepts=initial_accepts+1;
        if ls.increase && t<1 && norm(Fz0_new)>0
            increase_tests=increase_tests+1;
            t_hat=min(1,t/kappa);
            z0_hat=y0+t_hat*dk; Fz0_hat=PR.f(z0_hat,1); NF=NF+1; ls_nf=ls_nf+1;
            if line_test(PR,z0_hat,Fz0_hat,Fy0,NormFy0,dk,Normdk2,t_hat,para,ls,epsilon_C)
                t=t_hat; z0_new=z0_hat; Fz0_new=Fz0_hat;
                increase_accepts=increase_accepts+1;
            end
        end
    else
        for ls_it=1:ls_max
            tn=kappa*t;
            if tn<=0 || tn==t || isequal(z0_new,y0), break; end
            backtracks=backtracks+1;
            t=tn; z0_new=y0+t*dk; Fz0_new=PR.f(z0_new,1); NF=NF+1; ls_nf=ls_nf+1;
            ls_ok=line_test(PR,z0_new,Fz0_new,Fy0,NormFy0,dk,Normdk2,t,para,ls,epsilon_C);
            if ls_ok, break; end
        end
    end
    if ~ls_ok || isequal(z0_new,y0)
        L1=-1; failure_reason='line_search_failure'; break;
    end
    t_prev=t;
    accepted_steps(k)=t;

    NormFz0 = norm(Fz0_new);
    if NormFz0==0
        z1=y0;                   % F(z_k)=0 -> xi_k=0, x_{k+1}=P_C(varpi_k)
    elseif NormFz0>1e150
        % ||F(z)||^2 would overflow; use the normalised form.
        uz=Fz0_new/NormFz0;
        z1=y0-rho*dot(uz,y0-z0_new)*uz;
    else
        proj_coef=rho*dot(Fz0_new,y0-z0_new)/(NormFz0^2);
        z1=y0-proj_coef*Fz0_new;
    end
    x1  = PR.f(z1, 2);              % P_C
    Fx1 = PR.f(x1, 1);
    NF  = NF + 1;
    projection_nf=projection_nf+1;
    NormFx1 = norm(Fx1);
    if record_history
        hm(k)=signal_mse(x1,PR); ht(k)=toc(timer_id);
    end
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
assert(NF==1+inertia_nf+ls_nf+projection_nf,'Inconsistent NF accounting.');
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
    'actual_time',actual_time,'loop_iterations',k,'feasibility',feasibility,'x',x_return, ...
    'ls_nf',ls_nf,'backtracks',backtracks,'increase_tests',increase_tests, ...
    'increase_accepts',increase_accepts,'initial_accepts',initial_accepts, ...
    'inertia_nf',inertia_nf,'projection_nf',projection_nf, ...
    'accepted_steps',accepted_steps(isfinite(accepted_steps)));
if record_history
    hm(k)=signal_mse(x_return,PR); ht(k)=actual_time;
    info.mse=hm(1:k); info.time=ht(1:k);
end
if L1==1
    Itr=k; info.status='converged';
else
    NF=NaN; Itr=NaN; Tcpu=NaN; NormF=NaN;
end
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function ok=line_test(PR,z,Fz,Fw,nfw,d,nd2,t,p,ls,epsilon_C)
% Common exact-root stopping convention; see README.md.
ok=false;
if ~isreal(Fz) || any(~isfinite(Fz)), return; end
nfz=norm(Fz); if ~isfinite(nfz), return; end
% nfw is cached once at the inertial point and reused for every trial.
switch ls.criterion
    case 'Proposed'
        factor=max(p.mu,norm(Fz-Fw)/(p.a+nfz+nfw));
        rhs=ls.sigma*t*factor*nd2;
    case 'Yin2021'
        rhs=ls.sigma*t*max(ls.lambda,min(ls.nu,nfz))*nd2;
    case 'Liu'
        rhs=ls.sigma*t*nfz*min(ls.liu_mu,nd2);
    case 'Jiang'
        ndl=sqrt(nd2)^ls.j_l;
        fwl=nfw^(1-ls.j_l);
        etaJ=fwl/(1+ls.j_a*max(nfz,fwl*ndl));
        rhs=ls.sigma*t*etaJ*nfz*ndl;
end
lhs=-Fz'*d;
ok=isreal(lhs) && isfinite(lhs) && isreal(rhs) && isfinite(rhs) && lhs>0 && lhs>=rhs;
end

function value=signal_mse(x,PR)
value=mean((x(1:PR.signal_n)-x(PR.signal_n+1:end)-PR.true_x).^2);
end




