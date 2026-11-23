function res=huber_irls_same_model(inst,lam,target,max_iter,record_trace)
% Basic Huber IRLS for exactly min_(xi,e) .5||y-Gxi-e||^2+lam||e||_1.
% No Newton steps, ridge term, reweighted penalty, or change of loss.
if nargin<4,max_iter=2000;end
if nargin<5,record_trace=false;end
total_clock=tic;
G=inst.G;Gls=sparse(G);y=inst.y(:);m=inst.m;
assert(lam>0 && target>0 && max_iter>=1 && max_iter==fix(max_iter));
xi0=Gls\y;xi=xi0;Py=inst.proj(y);
if record_trace,history=nan(max_iter+1,3);else,history=zeros(0,3);end
setup_time=toc(total_clock);solver_clock=tic;
status='iteration_limit';it=0;checks=0;weighted_solves=0;monotonicity_violations=0;
first_weights=[];first_update=[];last_cost=Inf;
for k=0:max_iter
    r=y-Gls*xi;absr=abs(r);t=min(absr,lam);
    cost=sum(.5*t.^2+lam*(absr-t));
    if isfinite(last_cost) && cost>last_cost+64*eps(max(1,abs(last_cost)))
        monotonicity_violations=monotonicity_violations+1;
        status='objective_increase';break;
    end
    last_cost=cost;
    e=sign(r).*max(absr-lam,0);s=[max(e,0);max(-e,0)];
    Pe=inst.proj(e);
    fv=[min(s(1:m),Pe+lam-Py);min(s(m+1:end),-Pe+lam+Py)];
    fnorm=norm(fv);checks=checks+1;
    if record_trace,history(k+1,:)=[k,cost,fnorm];end
    if ~all(isfinite([xi;e;fv])) || ~isfinite(cost)
        status='nonfinite_irls_state';break;
    end
    if fnorm<=target,status='converged';break;end
    if k==max_iter,break;end
    weights=ones(m,1);out=absr>lam;weights(out)=lam./absr(out);
    sw=sqrt(weights);A=spdiags(sw,0,m,m)*Gls;b=y.*sw;
    next=A\b;weighted_solves=weighted_solves+1;
    if k==0,first_weights=weights;first_update=next;end
    if ~all(isfinite(next)),status='nonfinite_weighted_LS';break;end
    if isequal(next,xi),status='stagnation';break;end
    xi=next;it=it+1;
end
solver_wall=toc(solver_clock);recovery_clock=tic;
% Recompute from the actual last iterate even on a failure path.
r=y-Gls*xi;e=sign(r).*max(abs(r)-lam,0);s=[max(e,0);max(-e,0)];
xi_recovered=Gls\(y-e);fit=inst.proj(y-e);
objective=.5*norm(fit)^2+lam*norm(e,1);
split_objective=.5*norm(fit)^2+lam*sum(s);
t=min(abs(r),lam);huber_objective=sum(.5*t.^2+lam*(abs(r)-t));
relative_error=norm(xi_recovered-inst.true_xi)/norm(inst.true_xi);
recovery_time=toc(recovery_clock);verify_clock=tic;
% Operator form: P*(y-e) = Py - P*e, so no m-by-m matrix is ever formed.
Pe=inst.proj(e);
g=[lam*ones(m,1)-(Py-Pe);lam*ones(m,1)+(Py-Pe)];
natural=min(s,g);residual=norm(natural);feasibility=norm(min(s,0));
finite=all(isfinite([s;xi_recovered;natural])) && isfinite(objective);
success=strcmp(status,'converged') && finite && residual<=target && feasibility<=target;
if strcmp(status,'converged') && ~success,status='external_validation_failed';end
verification_time=toc(verify_clock);
if record_trace,history=history(1:checks,:);end
res=struct('success',success,'Status',status,'status',status,'Itr',it,'TotalItr',it, ...
    'NF',NaN,'actual_NF',NaN,'AttemptCount',1,'FirstAttemptSuccess',success, ...
    's',s,'e',e,'xi',xi_recovered,'xi_iterate',xi,'xi_initial',xi0, ...
    'RecoveryStateDifference',norm(xi_recovered-xi), ...
    'obj',objective,'split_obj',split_objective,'huber_objective',huber_objective, ...
    'relerr',relative_error,'residual',residual,'feasibility',feasibility, ...
    'native_stop_residual',fnorm,'WeightedLSSolves',weighted_solves,'CNEChecks',checks, ...
    'monotonicity_violations',monotonicity_violations,'history',history, ...
    'first_weights',first_weights,'first_update',first_update, ...
    'SetupTime',setup_time,'SolverWallTime',solver_wall,'RecoveryTime',recovery_time, ...
    'VerificationTime',verification_time,'TotalTime',toc(total_clock), ...
    'equation_tolerance',target,'feasibility_tolerance',target,'LSBackend','sparse-QR');
end
