function res=pmu_qp_reference(inst,lam,target,reference_options)
%PMU_QP_REFERENCE Independent convex-QP numerical accuracy reference.
% This is NOT a same-iteration-budget timing baseline or an exact certificate.
%
% min_(xi,e,r,t) .5*r'*r + lam*sum(t)
% s.t. G*xi + e + r = y,  e-t <= 0,  -e-t <= 0.
% Positive lam makes t=abs(e) at an optimum. The formulation is exactly
% min_(xi,e) .5*norm(y-G*xi-e)^2 + lam*norm(e,1); its Hessian is sparse PSD.
% No dense projector, Gram matrix, inverse, ridge term or altered loss is used.
%
% Native QP tolerances are distinct from the ORIGINAL equation tolerance.
% A positive exit flag alone never sets success=true. Returned points are
% checked using recovered-state residuals, QP KKT quantities and a separately
% projected/scaled numerical dual candidate. Floating-point checks do not
% establish an interval-arithmetic or exact mathematical certificate.
%
% Optional reference_options fields NativeTolerance and StepTolerance only
% control quadprog's internal stopping tests. Defaults preserve the original
% reference settings. Tighter native tolerances never relax or replace the
% independent target on the ORIGINAL CNE, and do not change the QP model.

total_clock=tic;
if nargin<4||isempty(reference_options),reference_options=struct();end
assert(isstruct(reference_options)&&isscalar(reference_options), ...
    'PMU:ReferenceOptions','reference_options must be a scalar struct.');
allowed_options={'NativeTolerance','StepTolerance'};
assert(all(ismember(fieldnames(reference_options),allowed_options)), ...
    'PMU:ReferenceOptions','Supported fields are NativeTolerance and StepTolerance only.');
if exist('quadprog','file')==0 || exist('optimoptions','file')==0 ...
        || ~license('test','Optimization_Toolbox')
    error('PMU:MissingOptimizationToolbox', ...
        ['pmu_qp_reference requires an installed and licensed MATLAB ' ...
         'Optimization Toolbox (quadprog). No substitute solver was run.']);
end
assert(isscalar(lam)&&isreal(lam)&&isfinite(lam)&&lam>0,'Invalid lambda.');
assert(isscalar(target)&&isreal(target)&&isfinite(target)&&target>0,'Invalid target.');
assert(all(isfield(inst,{'G','m','y','recover','proj'})), ...
    'The prepared PMU instance is missing required fields.');
assert(isa(inst.recover,'function_handle')&&isa(inst.proj,'function_handle'));
H=sparse(inst.G);y=inst.y(:);m=inst.m;nstate=size(H,2);
assert(isreal(H)&&isreal(y)&&size(H,1)==m&&numel(y)==m ...
    &&m>=nstate&&nstate>=1&&all(isfinite(nonzeros(H)))&&all(isfinite(y)), ...
    'Invalid real PMU measurement data.');

% x=[xi;e;r;t]. The identity in the equality matrix gives independent rows.
I=speye(m);Z=sparse(m,m);S=sparse(m,nstate);
nvar=nstate+3*m;
Q=spdiags([zeros(nstate+m,1);ones(m,1);zeros(m,1)],0,nvar,nvar);
linear=[zeros(nstate+2*m,1);lam*ones(m,1)];
Aeq=[H,I,I,Z];Aineq=[S,I,Z,-I;S,-I,Z,-I];bineq=zeros(2*m,1);

% These are solver-specific tolerances, not a claim that the original F norm
% is equal to a QP first-order or constraint metric.
native_tol=min(1e-10,max(1e-12,target*1e-4));
step_tol=1e-12;
if isfield(reference_options,'NativeTolerance')
    v=reference_options.NativeTolerance;
    assert(isnumeric(v)&&isscalar(v)&&isreal(v)&&isfinite(v)&&v>0, ...
        'PMU:ReferenceTolerance','NativeTolerance must be a positive finite real scalar.');
    native_tol=double(v);
end
if isfield(reference_options,'StepTolerance')
    v=reference_options.StepTolerance;
    assert(isnumeric(v)&&isscalar(v)&&isreal(v)&&isfinite(v)&&v>0, ...
        'PMU:ReferenceTolerance','StepTolerance must be a positive finite real scalar.');
    step_tol=double(v);
end
check_tol=max(target,100*native_tol);
agreement_tol=max(1e-10,target*1e-3);
max_iterations=1000;
options=optimoptions('quadprog','Algorithm','interior-point-convex', ...
    'Display','off','OptimalityTolerance',native_tol, ...
    'ConstraintTolerance',native_tol,'StepTolerance',step_tol, ...
    'MaxIterations',max_iterations,'LinearSolver','sparse');
setup_seconds=toc(total_clock);

res=struct('success',false,'status','quadprog_not_run', ...
    'ReferenceKind','convex-QP numerical accuracy reference', ...
    'ExactCertificate',false,'ReferenceOnly',true, ...
    'SameBudgetTimingBaseline',false, ...
    'TimingInterpretation','Reference construction and checks; not a same-budget speed comparison.', ...
    'exitflag',NaN,'output',struct(),'QPMultipliers',struct(), ...
    'Itr',NaN,'ActualItr',NaN,'NF',NaN,'ActualNF',NaN, ...
    'target',target,'OriginalTolerance',target,'feasibility_tolerance',target, ...
    'NativeOptimalityTolerance',native_tol,'NativeConstraintTolerance',native_tol, ...
    'NativeStepTolerance',step_tol,'RequestedReferenceOptions',reference_options, ...
    'NumericCheckTolerance',check_tol,'OriginalOracleAgreementTolerance',agreement_tol, ...
    'ReferenceMaxIterations',max_iterations,'SetupTime',setup_seconds, ...
    'SolverWallTime',NaN,'RecoveryTime',0,'VerificationTime',0,'TotalTime',NaN, ...
    'Residual',Inf,'residual',Inf,'feasibility',Inf,'Objective',NaN,'obj',NaN, ...
    'relerr',NaN,'s',[],'e',[],'xi',[],'xi_qp',[],'r_qp',[],'t_qp',[], ...
    'QPSolution',[],'QPReportedObjective',NaN,'QPRecomputedObjective',NaN, ...
    'QPObjectiveAgreement',NaN,'NativeRecoveredObjectiveDifference',NaN, ...
    'OriginalOracleAgreement',NaN,'OriginalTargetReached',false, ...
    'NumericKKTVerified',false,'NumericDualityVerified',false, ...
    'NumericalReferenceVerified',false,'KKT',struct(),'Duality',struct());

solver_clock=tic;
try
    [x,fval,flag,output,multipliers]=quadprog(Q,linear,Aineq,bineq,Aeq,y, ...
        [],[],[],options);
catch exception
    res.SolverWallTime=toc(solver_clock);
    res.status='quadprog_exception';
    res.output=struct('message',exception.message,'ExceptionIdentifier',exception.identifier);
    res.TotalTime=toc(total_clock);
    return;
end
res.SolverWallTime=toc(solver_clock);
res.exitflag=flag;res.output=output;res.QPMultipliers=multipliers;
res.QPReportedObjective=fval;
if isfield(output,'iterations'),res.Itr=output.iterations;res.ActualItr=output.iterations;end
if isempty(x)||numel(x)~=nvar||~isreal(x)||any(~isfinite(x))
    res.status='quadprog_no_finite_reference_point';res.TotalTime=toc(total_clock);return;
end
x=x(:);res.QPSolution=x;
xi_qp=x(1:nstate);e=x(nstate+(1:m));r_qp=x(nstate+m+(1:m));
t_qp=x(nstate+2*m+(1:m));

recovery_clock=tic;
s=[max(e,0);max(-e,0)];xi=inst.recover(y-e);
xi=xi(:);assert(numel(xi)==nstate,'Prepared recovery returned a wrong state size.');
fit=y-e-H*xi;
objective=.5*dot(fit,fit)+lam*norm(e,1);
res.RecoveryTime=toc(recovery_clock);

verify_clock=tic;
% Independent original CNE check through a recovered LS fit. In exact
% arithmetic fit=P(y-e), so the two blocks below equal min(s,K*s+c).
original=[min(s(1:m),lam-fit);min(s(m+1:end),lam+fit)];
residual=norm(original);feasibility=norm(min(s,0));
Pe=inst.proj(e);
if isfield(inst,'Py'),Py=inst.Py;else,Py=inst.proj(y);end
Pe=Pe(:);Py=Py(:);
assert(numel(Pe)==m&&numel(Py)==m,'Prepared projector returned a wrong size.');
oracle=[min(s(1:m),Pe+lam-Py);min(s(m+1:end),-Pe+lam+Py)];
agreement=norm(original-oracle);
qp_objective=.5*dot(r_qp,r_qp)+lam*sum(t_qp);

% Independent KKT residuals for exactly the lifted QP supplied to quadprog.
kkt=struct('Available',false,'EqualityInf',Inf,'EqualityScaled',Inf, ...
    'InequalityInf',Inf,'InequalityScaled',Inf,'DualInequalityInf',Inf, ...
    'DualInequalityScaled',Inf,'StationarityInf',Inf,'StationarityScaled',Inf, ...
    'ComplementarityInf',Inf,'ComplementarityScaled',Inf,'ResidualScaled',Inf);
have_duals=isstruct(multipliers)&&isfield(multipliers,'eqlin') ...
    &&isfield(multipliers,'ineqlin')&&numel(multipliers.eqlin)==m ...
    &&numel(multipliers.ineqlin)==2*m;
if have_duals
    nu=multipliers.eqlin(:);mu=multipliers.ineqlin(:);
    have_duals=isreal(nu)&&isreal(mu)&&all(isfinite([nu;mu]));
end
if have_duals
    eq=Aeq*x-y;ineq=Aineq*x-bineq;slack=-ineq;
    gradient=Q*x+linear;di=Aineq'*mu;de=Aeq'*nu;
    stationarity=gradient+di+de;
    eq_scale=max([1,norm(y,inf),norm(Aeq*x,inf)]);
    ineq_scale=max(1,norm(Aineq*x,inf));
    dual_scale=max(1,norm(mu,inf));
    stat_scale=max([1,norm(gradient,inf),norm(di,inf),norm(de,inf)]);
    comp_scale=max(1,abs(qp_objective));
    kkt.Available=true;kkt.EqualityInf=norm(eq,inf);
    kkt.EqualityScaled=kkt.EqualityInf/eq_scale;
    kkt.InequalityInf=max([0;ineq]);kkt.InequalityScaled=kkt.InequalityInf/ineq_scale;
    kkt.DualInequalityInf=max([0;-mu]);
    kkt.DualInequalityScaled=kkt.DualInequalityInf/dual_scale;
    kkt.StationarityInf=norm(stationarity,inf);
    kkt.StationarityScaled=kkt.StationarityInf/stat_scale;
    kkt.ComplementarityInf=norm(mu.*slack,inf);
    kkt.ComplementarityScaled=kkt.ComplementarityInf/comp_scale;
    kkt.ResidualScaled=max([kkt.EqualityScaled,kkt.InequalityScaled, ...
        kkt.DualInequalityScaled,kkt.StationarityScaled,kkt.ComplementarityScaled]);
    dual_raw=-nu;dual_source='negative QP equality multiplier';
else
    dual_raw=fit;dual_source='recovered LS residual (QP multipliers unavailable)';
end

% Original-model dual: max y'*v-.5*v'*v, G'*v=0, norm(v,inf)<=lam.
% Projection enforces range(H)^perp in exact arithmetic; positive scalar
% contraction preserves it and enforces the box. Numerical residuals remain
% explicit: this is not an exact floating-point feasibility certificate.
dual_projected=inst.proj(dual_raw);dual_projected=dual_projected(:);
assert(numel(dual_projected)==m,'Prepared projector returned a wrong dual size.');
dual_scale=1;dual_max=norm(dual_projected,inf);
if dual_max>lam,dual_scale=(lam/dual_max)*(1-16*eps);end
dual=dual_scale*dual_projected;
dual_eq=norm(H'*dual,inf);dual_box=max(0,norm(dual,inf)-lam);
dual_eq_scale=max(1,norm(H,1)*norm(dual,inf));
dual_value=dot(y,dual)-.5*dot(dual,dual);
gap=objective-dual_value;gap_scale=max([1,abs(objective),abs(dual_value)]);
scaled_gap=gap/gap_scale; % Scaled gap, not a relative gap when objectives <1.
dual_check=max(dual_eq/dual_eq_scale,dual_box/max(1,lam));
duality=struct('Source',dual_source,'ProjectedCandidate',true, ...
    'Scale',dual_scale,'Vector',dual,'EqualityInf',dual_eq, ...
    'EqualityScaled',dual_eq/dual_eq_scale,'BoxViolationInf',dual_box, ...
    'FeasibilityScaled',dual_check,'PrimalObjective',objective, ...
    'DualObjective',dual_value,'Gap',gap,'GapScale',gap_scale, ...
    'ScaledGap',scaled_gap,'Tolerance',check_tol, ...
    'ExactFeasibilityCertified',false);

finite=all(isfinite([s;xi;original;oracle;dual])) ...
    &&all(isfinite([objective,qp_objective,dual_value,gap,agreement]));
objective_agreement=Inf;
if isscalar(fval)&&isreal(fval)&&isfinite(fval)
    objective_agreement=abs(fval-qp_objective)/max([1,abs(fval),abs(qp_objective)]);
end
original_pass=finite&&residual<=target&&feasibility<=target&&agreement<=agreement_tol;
kkt_pass=kkt.Available&&isfinite(kkt.ResidualScaled)&&kkt.ResidualScaled<=check_tol;
duality_pass=isfinite(dual_check)&&isfinite(scaled_gap)&&dual_check<=check_tol ...
    &&abs(scaled_gap)<=check_tol;
verified=flag>0&&original_pass&&kkt_pass&&duality_pass&&objective_agreement<=check_tol;

res.s=s;res.e=e;res.xi=xi;res.xi_qp=xi_qp;res.r_qp=r_qp;res.t_qp=t_qp;
res.obj=objective;res.Objective=objective;res.residual=residual;res.Residual=residual;
res.feasibility=feasibility;res.OriginalOracleAgreement=agreement;
res.QPRecomputedObjective=qp_objective;res.QPObjectiveAgreement=objective_agreement;
res.NativeRecoveredObjectiveDifference=qp_objective-objective;
if isfield(inst,'true_xi')
    truth=inst.true_xi(:);
    if numel(truth)==nstate&&isreal(truth)&&all(isfinite(truth))&&norm(truth)>0
        res.relerr=norm(xi-truth)/norm(truth);
    end
end
res.KKT=kkt;res.Duality=duality;res.OriginalTargetReached=original_pass;
res.NumericKKTVerified=kkt_pass;res.NumericDualityVerified=duality_pass;
res.NumericalReferenceVerified=verified;res.success=verified;
if verified,res.status='numerical_reference_verified';
elseif flag<=0,res.status='quadprog_exit_not_success';
elseif ~original_pass,res.status='original_equation_checks_not_passed';
else,res.status='numerical_optimality_checks_not_passed';end
res.VerificationTime=toc(verify_clock);res.TotalTime=toc(total_clock);
end
