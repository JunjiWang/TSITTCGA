function res=pmu_run_cgpm(inst,solver,p,lam,initial)
% Common adapter: original F, native method, independent endpoint and NF check.
if nargin<5,initial='raw';end
t0=tic;m=inst.m;Py=inst.Py;
if strcmp(initial,'raw'),e0=Py;
elseif strcmp(initial,'soft'),e0=sign(Py).*max(abs(Py)-lam,0);
else,error('PMU:Initial','Use raw or soft initialization.');end
x0=[max(e0,0);max(-e0,0)];calls=0;projections=0;
PR=struct('name','PMU','n',2*m,'f',@oracle);
setup=toc(t0);clock=tic;
[tc,nf,it,reported,info]=solver(PR,x0,p);
wall=toc(clock);s=info.x(:);assert(numel(s)==2*m,'Invalid endpoint dimension.');
clock=tic;e=s(1:m)-s(m+1:end);xi=inst.recover(inst.y-e);
% Independent LS-based residual, rather than the oracle's difference P*e-P*y.
fit=inst.y-e-inst.G*xi;
fv=[min(s(1:m),lam-fit);min(s(m+1:end),lam+fit)];
residual=norm(fv);feasibility=norm(min(s,0));
obj=.5*dot(fit,fit)+lam*norm(e,1);
finite=all(isfinite([s;xi;fv]))&&isfinite(obj);
count_ok=calls==info.actual_NF;
if isfinite(nf),count_ok=count_ok&&nf==calls;end
% NormF is NaN on some native failure paths; do not reject a failure record.
report_ok=~info.success || (isfinite(reported)&&abs(reported-residual)<=1e-9);
reached=finite&&residual<=p.epsilon&&feasibility<=p.epsilon_C;
success=logical(info.success)&&reached&&count_ok&&report_ok;
status=info.status;
if info.success&&~success,status='external_validation_failed';end
assert(count_ok,'PMU:NF','Native and counted NF differ.');
assert(report_ok,'PMU:Residual','Reported successful residual differs.');
validation=toc(clock);
floor_accepts=0;
if isfield(info,'paper_floor_accepts'),floor_accepts=info.paper_floor_accepts;end
res=struct('success',success,'status',status,'native_success',logical(info.success), ...
    'criterion_reached',reached,'Itr',it,'ActualItr',info.loop_iterations, ...
    'NF',nf,'ActualNF',info.actual_NF,'NFCounted',calls,'ProjectionCount',projections, ...
    'count_ok',count_ok,'report_ok',report_ok,'residual',residual,'feasibility',feasibility, ...
    'obj',obj,'split_obj',.5*dot(fit,fit)+lam*sum(s),'relerr',norm(xi-inst.true_xi)/norm(inst.true_xi), ...
    'SetupTime',setup,'SolverWallTime',wall,'ValidationRecoveryTime',validation, ...
    'TotalTime',toc(t0),'NativeTime',tc,'s',s,'e',e,'xi',xi,'x0',x0, ...
    'WeightedLSSolves',NaN,'CNEChecks',NaN,'PaperFloorAccepts',floor_accepts,'Initial',initial,'info',info);
    function v=oracle(x,mode)
        if mode==1,calls=calls+1;elseif mode==2,projections=projections+1;end
        v=pmu_map(x,mode,inst,lam);
    end
end
