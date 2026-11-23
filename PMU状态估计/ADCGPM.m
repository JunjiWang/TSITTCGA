function [Tcpu,NF,Itr,NormF,info] = ADCGPM(PR,x0,para)
% Jiang, Luo and Jian, JOTA 210:64 (2026), Algorithm 1.
% DOI 10.1007/s10957-026-03090-3. Defaults reproduce Section 5.1.2:
% p_k=d_{k-1}, q_k=F_k, with no inertial extrapolation.
% Eq. (19) uses alpha=shrink^i, i=0,1,..., starting from alpha=1.
% Undefined denominators are reported, never replaced by a new direction.
if nargin<3 || isempty(para), para=struct; end
clock_id=tic; % Include parameter checks and initial feasibility checks, as in the other solvers.
p=struct('ell',0.5,'sigma',0.01,'shrink',0.55,'rho1',10,'rho2',50, ...
    'mu',0.4,'zeta',50,'xi',0.01,'a',2,'gamma',1.9, ...
    'Itr_max',2000,'eps_F',1e-6, ...
    'ls_max',2000,'quiet',true);
fields=fieldnames(p);
for j=1:numel(fields)
    if isfield(para,fields{j}) && ~isempty(para.(fields{j}))
        p.(fields{j})=para.(fields{j});
    end
end
x=x0(:);
Fx=PR.f(x,1); NF=1; fn=norm(Fx);
status='iteration_limit'; success=false;
x_prev=[]; F_prev=[]; d_prev=[];
updates=0; max_backtracks=0;
for k=0:p.Itr_max-1
    if ~isreal(Fx) || any(~isfinite(Fx)) || ~isfinite(fn)
        status='nonfinite_function_value'; break;
    end
    if fn<=p.eps_F
        success=true; status='converged'; break;
    end
    if k==0
        d=-Fx;
    else
        y=Fx-F_prev;
        if norm(y)==0
            d=-Fx;
        else
            s=x-x_prev;
            sy=s'*y; dy=d_prev'*y;
            if ~isfinite(sy) || ~isfinite(dy)
                status='nonfinite_secant_denominator'; break;
            elseif sy==0
                status='undefined_secant_curvature'; break;
            elseif dy==0
                status='undefined_conjugate_denominator'; break;
            end
            % Eq. (8), with the original unmodified secant pair.
            beta=2*(Fx'*y)/dy-(1+(y'*y)/sy)*(Fx'*s)/dy;
            if ~isreal(beta) || ~isfinite(beta)
                status='nonfinite_conjugate_parameter'; break;
            end
            % Eq. (10), q_k=F_k in Section 5.1.2.
            q=Fx;
            clip=abs(q'*Fx)/(norm(q)*norm(d_prev));
            if ~isfinite(clip)
                status='undefined_truncation_denominator'; break;
            end
            beta_star=min(p.rho1*clip,max(beta,-p.rho2*clip));
            cross=beta_star*(Fx'*d_prev);
            if cross>-p.zeta*fn^2 && cross<p.mu*fn^2
                d=-Fx+beta_star*d_prev;
            else
                % Eq. (11), p_k=d_{k-1} in Section 5.1.2.
                pk=d_prev; pnorm2=pk'*pk;
                if pnorm2==0 || ~isfinite(pnorm2)
                    status='undefined_restart_denominator'; break;
                end
                d=-Fx+p.xi*((Fx'*pk)/pnorm2)*pk;
            end
        end
    end
    nd=norm(d);
    if ~isreal(d) || any(~isfinite(d)) || ~isfinite(nd)
        status='nonfinite_direction'; break;
    end
    % Eqs. (19)--(20), fixed geometric trials 1,shrink,shrink^2,... .
    fixed=fn^(1-p.ell); nd_l=nd^p.ell;
    if ~isfinite(fixed) || ~isfinite(nd_l)
        status='nonfinite_line_search_scaling'; break;
    end
    accepted=false;
    for i=0:p.ls_max
        alpha=p.shrink^i;
        z=x+alpha*d;
        Fz=PR.f(z,1); NF=NF+1; fz=norm(Fz);
        eta=fixed/(1+p.a*max(fz,fixed*nd_l));
        lhs=-Fz'*d;
        rhs=p.sigma*alpha*eta*fz*nd_l;
        accepted=isreal(Fz) && all(isfinite(Fz)) && isfinite(fz) ...
            && isreal(lhs) && isfinite(lhs) && isreal(rhs) && isfinite(rhs) && lhs>=rhs;
        if accepted, break; end
        if alpha<=0 || isequal(z,x), break; end
    end
    max_backtracks=max(max_backtracks,i);
    if ~accepted || isequal(z,x)
        status='line_search_failure'; break;
    end
    % Eq. (22): F(z_k)=0 gives lambda_k=0, so no division by zero here.
    if fz == 0
        lambda = 0;
    else
        lambda = (Fz'*(x-z))/fz^2;
    end
    x_next=PR.f(x-p.gamma*lambda*Fz,2);  % Eq. (21).
    if ~isreal(x_next) || any(~isfinite(x_next))
        status='nonfinite_projection'; break;
    end
    x_prev=x; F_prev=Fx; d_prev=d;
    x=x_next; updates=updates+1;
    Fx=PR.f(x,1); NF=NF+1; fn=norm(Fx);
    if ~isreal(Fx) || any(~isfinite(Fx)) || ~isfinite(fn)
        status='nonfinite_projected_function_value'; break;
    end
    if fn<=p.eps_F
        success=true; status='converged'; break;
    end
end
Tcpu=toc(clock_id);
actual_NF=NF; actual_time=Tcpu;
feasibility=norm(x-PR.f(x,2));
info=struct('success',success,'status',status,'actual_NF',actual_NF, ...
    'actual_time',actual_time,'loop_iterations',k+1,'x',x, ...
    'feasibility',feasibility,'accepted_updates',updates, ...
    'max_backtracks',max_backtracks);
if success
    Itr=k+1; NormF=fn;
else
    Tcpu=NaN; NF=NaN; Itr=NaN; NormF=NaN;
end
if ~p.quiet
    fprintf('ADCGPM: %s, NF=%d, updates=%d, ||F||=%.4e\n', ...
        status,actual_NF,updates,fn);
end
end

