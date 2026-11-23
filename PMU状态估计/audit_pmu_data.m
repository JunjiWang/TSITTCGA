function report=audit_pmu_data()
%AUDIT_PMU_DATA Check frozen PMU data and the split-CNE implementation only.
% No estimation algorithm or optimization solver is run. Identical G/Q
% pairs share a matrix check; every observation/noise/oracle is still checked.
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');
if ~isfolder(out),mkdir(out);end
oldrng=rng;oldpath=path;cleanup=onCleanup(@()restore_session(oldrng,oldpath)); %#ok<NASGU>
addpath(here,'-begin');
actual=which('pmu_map');
assert(strcmpi(strrep(actual,'/','\'),strrep(fullfile(here,'pmu_map.m'),'/','\')), ...
    'PMUDataAudit:ShadowedMap','A different pmu_map shadows this package.');
lam=.005;cache=struct([]);rows=struct([]);file_rows=struct([]);
systems={'IEEE14','IEEE118','IEEE300'};buses=[14 118 300];
for grid=1:3
    file=fullfile(here,'data',sprintf('input_ieee%d.mat',buses(grid)));
    try
        d=load(file);
        if isfield(d,'instance_records'),instances=d.instance_records;
        elseif isfield(d,'instances'),instances=d.instances;
        else,error('PMUDataAudit:MissingInstances','No instance_records or instances variable.');end
        need(iscell(instances)&&numel(instances)==45,'Expected 45 cell-array instances per network.');
        for key={'lambda','lam'}
            if isfield(d,key{1})&&isnumeric(d.(key{1}))&&isscalar(d.(key{1}))
                need(d.(key{1})==lam,'Stored file-level regularization differs from .005.');
            end
        end
        circuit=[];
        if isfield(d,'c')&&isstruct(d.c)&&isscalar(d.c)&&isfield(d.c,'branch'),circuit=d.c;end
        if isfield(d,'network_check')&&isstruct(d.network_check)&&isfield(d.network_check,'passed')
            need(d.network_check.passed,'Stored physical-network verification did not pass.');
        end
        for k=1:numel(instances)
            row=empty_row(systems{grid},k);
            row.CircuitMetadataAvailable=~isempty(circuit)&&all(isfield(circuit.branch,{'f','t','Yff','Yft','Ytf','Ytt'}));
            try
                a=instances{k};
                required={'G','H','y','ytrue','true_xi','bad','noise','n','nq','m', ...
                    'sigma','seed','nbad','bad_indices','bus_selection','current_selection'};
                need(isstruct(a)&&isscalar(a)&&all(isfield(a,required)),'Missing pmu_build data fields.');
                n=a.n;m=a.m;G=a.G;H=a.H;Q=[];
                if isfield(a,'Q'),Q=a.Q;end
                row.StoredQ=~isempty(Q);
                need(n==buses(grid)&&m==2*a.nq&&isequal(size(G),[m 2*n]) ...
                    &&isequal(size(H),[a.nq n])&&m>=2*n,'Inconsistent real/complex measurement dimensions.');
                need(isreal(G)&&all(isfinite(nonzeros(G)))&&all(isfinite(nonzeros(H))) ...
                    &&all(isfinite([a.y(:);a.ytrue(:);a.true_xi(:);a.bad(:);a.noise(:)])), ...
                    'Nonfinite measurement matrix, physical state or observations.');
                need(numel(a.y)==m&&numel(a.ytrue)==m&&numel(a.bad)==m&&numel(a.noise)==m ...
                    &&numel(a.true_xi)==2*n,'Observation/state vector shape mismatch.');
                need(a.sigma==.001&&a.nbad==nnz(a.bad)&&a.nbad==numel(a.bad_indices) ...
                    &&numel(unique(a.bad_indices))==a.nbad,'Noise level or bad-data count mismatch.');
                row.StateDimension=2*n;row.MeasurementDimension=m;row.CNEDimension=2*m;
                row.Currents=numel(a.current_selection);row.BadDataCount=a.nbad;
                row.StoredLambda=isfield(a,'lambda');
                if row.StoredLambda,need(a.lambda==lam,'Stored instance lambda differs from .005.');end
                if isfield(a,'Currents'),need(a.Currents==row.Currents,'Current selection count mismatch.');end
                if isfield(a,'BadData'),need(a.BadData==a.nbad,'Stored BadData count differs.');end
                if isfield(a,'RandomSeed'),need(a.RandomSeed==a.seed,'Stored seed aliases differ.');end
                row.Seed=a.seed;row.NoiseIdentityError=norm(a.y(:)-a.ytrue(:)-a.noise(:)-a.bad(:),inf);
                need(row.NoiseIdentityError<=1e-12,'y=ytrue+noise+bad identity failed.');
                row.TruthIdentityError=norm(G*a.true_xi(:)-a.ytrue(:),inf);
                need(row.TruthIdentityError<=1e-10*max(1,norm(a.ytrue,inf)),'G*true_xi differs from ytrue.');
                Greal=[real(H),-imag(H);imag(H),real(H)];
                row.RealComplexMatrixError=norm(G-Greal,'fro')/max(1,norm(G,'fro'));
                V=a.true_xi(1:n)+1i*a.true_xi(n+1:end);hv=H*V;
                row.ComplexTruthError=norm([real(hv);imag(hv)]-a.ytrue(:),inf);
                need(row.RealComplexMatrixError<=1e-12 ...
                    &&row.ComplexTruthError<=1e-10*max(1,norm(a.ytrue,inf)), ...
                    'Complex PMU H and real G/state/measurement ordering differ.');
                buses_selected=sort(a.bus_selection(:));
                need(numel(unique(buses_selected))==numel(buses_selected) ...
                    &&all(buses_selected>=1&buses_selected<=n&buses_selected==fix(buses_selected)), ...
                    'Invalid selected voltage buses.');
                nv=numel(buses_selected);voltage=sparse((1:nv)',buses_selected,ones(nv,1),nv,n);
                row.VoltageRowsError=norm(H(1:nv,:)-voltage,'fro');
                need(row.VoltageRowsError==0,'Selected voltage rows are not the identity-row measurements.');
                row.ExactVoltageRankCertificate=nv==n;
                if row.ExactVoltageRankCertificate
                    need(nnz(G([1:n,a.nq+(1:n)],:)-speye(2*n))==0,'All-voltage rank certificate failed.');
                end
                rng(a.seed,'twister');bad=zeros(m,1);idx=[];
                if a.nbad>0
                    idx=randperm(m,a.nbad);
                    bad(idx)=(20+80*rand(a.nbad,1)).*sign(randn(a.nbad,1))*a.sigma;
                end
                noise=a.sigma*randn(m,1);
                row.SeedRegenerationExact=isequal(bad,a.bad(:))&&isequal(noise,a.noise(:)) ...
                    &&isequal(idx(:),a.bad_indices(:));
                need(row.SeedRegenerationExact,'Recorded pmu_build seed does not reproduce bad data/noise.');
                if row.StoredQ
                    need(isequal(size(Q),[m 2*n])&&isreal(Q)&&all(isfinite(Q(:))),'Stored economy Q is malformed.');
                end
                key=0;
                for j=1:numel(cache)
                    if isequal(G,cache(j).G)&&(~row.StoredQ||isequal(Q,cache(j).Q))
                        key=j;Q=cache(j).Q;break;
                    end
                end
                if key==0
                    if ~row.StoredQ,[Q,~]=qr(G,0);end
                    entry=check_matrix(G,Q);cache=push(cache,entry);key=numel(cache);
                end
                entry=cache(key);row.MatrixKey=key;
                for name={'QOrthogonalityError','QRangeError','RTriangularError','RReciprocalCondition', ...
                        'ProjectorIdempotenceError','ProjectorAdjointError','IndependentProjectorError'}
                    row.(name{1})=entry.(name{1});
                end
                row.NumericalRank=entry.Rank;
                if isfield(a,'rank'),need(a.rank==2*n,'Stored rank differs from verified full column rank.');end
                % Use the same projector backend as the solver, so this audit checks the operator the experiments actually run with (normal equations with a sparse Cholesky factor).  The explicit Q computed above is still checked for orthogonality and range.
                cc=pmu_parameters();[inst,~]=pmu_prepare_backend(a,cc.backend);
                py_reference=a.y(:)-entry.Gs*(entry.Factor\a.y(:));
                row.ProjectedObservationError=norm(inst.Py-py_reference)/max(1,norm(py_reference));
                need(row.ProjectedObservationError<=1e-9,'Projected observations disagree with independent sparse LS.');
                [row.OracleAgreementError,row.MonotonicityMinimum,row.LipschitzMaximumRatio, ...
                    row.OrthantProjectionError]=check_oracle(inst,entry,lam);
                if ~isempty(circuit)
                    [row.CircuitRowsChecked,row.CircuitRowsError]=check_circuit_rows(a,circuit);
                end
                row.Pass=true;
            catch exception
                row.Pass=false;row.ErrorIdentifier=exception.identifier;row.ErrorMessage=exception.message;
            end
            rows=push(rows,row);
        end
        file_rows=push(file_rows,struct('System',systems{grid},'Cases',numel(instances),'Loaded',true,'Error',''));
    catch exception
        row=empty_row(systems{grid},0);row.ErrorIdentifier=exception.identifier;row.ErrorMessage=exception.message;
        rows=push(rows,row);file_rows=push(file_rows,struct('System',systems{grid},'Cases',0,'Loaded',false,'Error',exception.message));
    end
end
CaseRows=struct2table(rows);valid=CaseRows.Case>0;
Summary=struct('CheckedLocal',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'ExpectedCases',135,'ObservedCases',nnz(valid),'PassedCases',nnz(valid&CaseRows.Pass), ...
    'FailedCases',nnz(valid&~CaseRows.Pass),'FileLoadFailures',nnz(~[file_rows.Loaded]), ...
    'FixedLambda',lam,'NoiseStdPerRealComponent',.001,'UniqueMeasurementMatrixChecks',numel(cache), ...
    'CasesWithStoredLambda',nnz(valid&CaseRows.StoredLambda), ...
    'CasesWithStoredQ',nnz(valid&CaseRows.StoredQ), ...
    'ComplexHChecks',nnz(valid&isfinite(CaseRows.RealComplexMatrixError)), ...
    'CircuitCoefficientRowChecks',nnz(valid&CaseRows.CircuitRowsChecked), ...
    'CircuitMetadataUnavailableCases',nnz(valid&~CaseRows.CircuitMetadataAvailable), ...
    'MaxQOrthogonalityError',finite_max(CaseRows.QOrthogonalityError), ...
    'MaxQRangeError',finite_max(CaseRows.QRangeError), ...
    'MaxIndependentProjectorError',finite_max(CaseRows.IndependentProjectorError), ...
    'MaxOracleAgreementError',finite_max(CaseRows.OracleAgreementError), ...
    'MinimumSampleMonotonicity',finite_min(CaseRows.MonotonicityMinimum), ...
    'MaximumSampleLipschitzRatio',finite_max(CaseRows.LipschitzMaximumRatio), ...
    'Scope','Frozen data, available circuit fields and implementation checks; no estimation solver executed', ...
    'TheoryBoundary','Finite samples of monotonicity/2-Lipschitz continuity validate implementation only; they are not a proof.', ...
    'UnavailableMetadata','Missing file-level circuit coefficients are reported as unchecked; H/G/complex truth are still checked.', ...
    'Files',file_rows,'Pass',nnz(valid)==135&&all(CaseRows.Pass)&&all([file_rows.Loaded]));
report=struct('Summary',Summary,'CaseRows',CaseRows,'Pass',Summary.Pass);
save(fullfile(out,'data_audit.mat'),'report','Summary','CaseRows','-v7');
writetable(CaseRows,fullfile(out,'data_audit.csv'));
short=struct('Summary',Summary,'FailedCaseRows',table2struct(CaseRows(~CaseRows.Pass,:)));
write_json(fullfile(out,'data_audit.json'),short);disp(Summary);
end

function entry=check_matrix(G,Q)
[m,n]=size(G);entry=struct('G',G,'Q',Q,'Gs',sparse(G),'Factor',[], ...
    'QOrthogonalityError',norm(Q'*Q-eye(n),'fro'), ...
    'QRangeError',norm(G-Q*(Q'*G),'fro')/max(1,norm(G,'fro')), ...
    'RTriangularError',NaN,'RReciprocalCondition',NaN,'Rank',0, ...
    'ProjectorIdempotenceError',NaN,'ProjectorAdjointError',NaN,'IndependentProjectorError',NaN);
need(entry.QOrthogonalityError<=1e-9&&entry.QRangeError<=1e-9,'Stored Q is not an orthonormal complete range basis.');
R=Q'*G;entry.RTriangularError=norm(tril(R,-1),'fro')/max(1,norm(R,'fro'));
entry.RReciprocalCondition=rcond(full(R));tol=max(m,n)*eps(max(1,norm(R,inf)));
entry.Rank=nnz(abs(diag(R))>tol);
need(entry.RTriangularError<=1e-9&&entry.Rank==n&&entry.RReciprocalCondition>1e-12, ...
    'Economy QR does not certify numerically full column rank.');
P=@(v)v-Q*(Q'*v);u=sin((1:m)'*sqrt(2));u=u/norm(u);
v=cos((1:m)'*sqrt(3));v=v/norm(v);w=sin((1:n)'*sqrt(5));w=w/norm(w);
probes=[u,v,G*w];projected=P(probes);
reference=probes-entry.Gs*(entry.Gs\probes);
entry.IndependentProjectorError=norm(projected-reference,'fro')/max(1,norm(reference,'fro'));
entry.ProjectorIdempotenceError=norm(P(projected)-projected,'fro')/max(1,norm(projected,'fro'));
entry.ProjectorAdjointError=abs(u'*P(v)-P(u)'*v);
need(entry.IndependentProjectorError<=1e-9&&entry.ProjectorIdempotenceError<=1e-9 ...
    &&entry.ProjectorAdjointError<=1e-9,'Projector failed independent LS/idempotence/self-adjoint checks.');
entry.Factor=decomposition(entry.Gs,'qr');
end
function [agreement,monotonicity,lip,projection_error]=check_oracle(inst,entry,lam)
m=inst.m;j=(1:2*m)';agreement=0;monotonicity=Inf;lip=0;projection_error=0;
for trial=1:3
    if trial==1,s=.02*sin(j*sqrt(2));t=.03*cos(j*sqrt(3));
    elseif trial==2,s=max(.04*sin(j*sqrt(5)),0);t=max(.03*cos(j*sqrt(7)),0);
    else,s=sin(j*.13);t=.5*cos(j*.17);end
    fs=pmu_map(s,1,inst,lam);ft=pmu_map(t,1,inst,lam);
    e=s(1:m)-s(m+1:end);r=inst.y(:)-e;
    fit=r-entry.Gs*(entry.Factor\r);
    independent=[min(s(1:m),lam-fit);min(s(m+1:end),lam+fit)];
    agreement=max(agreement,norm(fs-independent)/max(1,norm(independent)));
    ds=s-t;df=fs-ft;value=dot(df,ds);scale=max(1,norm(df)*norm(ds));
    monotonicity=min(monotonicity,value);ratio=norm(df)/norm(ds);lip=max(lip,ratio);
    need(all(isfinite([fs;ft]))&&value>=-1e-9*scale&&ratio<=2+1e-9, ...
        'Split-CNE implementation fails a monotonicity or 2-Lipschitz sample.');
    projected=pmu_map(s,2,inst,lam);
    projection_error=max(projection_error,norm(projected-max(s,0)));
end
need(agreement<=1e-9&&projection_error==0,'PMU oracle or orthant projection implementation differs.');
end
function [checked,error_value]=check_circuit_rows(a,c)
checked=false;error_value=NaN;E=c.branch;
keys={'f','t','Yff','Yft','Ytf','Ytt'};
if ~all(isfield(E,keys)),return;end
n=a.n;if isfield(c,'n'),need(c.n==n,'Circuit and instance bus counts differ.');end
nl=numel(E.f);need(numel(E.t)==nl,'Circuit branch endpoint lengths differ.');
sel=false(n,1);sel(a.bus_selection)=true;expected=sparse(0,n);
for i=find(sel)'
    row=sparse(1,i,1,1,n);expected=[expected;row]; %#ok<AGROW>
end
for k=1:nl
    f=E.f(k);t=E.t(k);
    if sel(f)&&ismember(2*k-1,a.current_selection)
        row=sparse([1 1],[f t],[E.Yff(k) E.Yft(k)],1,n);expected=[expected;row]; %#ok<AGROW>
    end
    if sel(t)&&ismember(2*k,a.current_selection)
        row=sparse([1 1],[f t],[E.Ytf(k) E.Ytt(k)],1,n);expected=[expected;row]; %#ok<AGROW>
    end
end
need(isequal(size(expected),size(a.H)),'Circuit branch-selection row count differs from H.');
error_value=norm(a.H-expected,'fro')/max(1,norm(a.H,'fro'));
need(error_value<=1e-12,'PMU H current rows differ from the stored complex branch coefficients.');
% When physical pi-model parameters exist, also check those coefficients.
if all(isfield(E,{'r','x','b','tap','shift','status'}))
    tap=E.tap;tap(tap==0)=1;tap=tap.*exp(1i*pi/180*E.shift);
    ys=E.status./(E.r+1i*E.x);bc=1i*E.status.*E.b/2;
    errors=[norm(E.Yff-(ys+bc)./abs(tap).^2,inf),norm(E.Yft+ys./conj(tap),inf), ...
        norm(E.Ytf+ys./tap,inf),norm(E.Ytt-ys-bc,inf)];
    scale=max([1;abs(E.Yff(:));abs(E.Yft(:));abs(E.Ytf(:));abs(E.Ytt(:))]);
    need(all(isfinite(errors))&&max(errors)<=1e-10*scale,'Stored pi-model current coefficients are inconsistent.');
end
checked=true;
end
function r=empty_row(system,k)
r=struct('System',system,'Case',k,'StateDimension',NaN,'MeasurementDimension',NaN,'CNEDimension',NaN, ...
    'Currents',NaN,'BadDataCount',NaN,'Seed',NaN,'StoredLambda',false,'StoredQ',false, ...
    'NoiseIdentityError',NaN,'TruthIdentityError',NaN,'RealComplexMatrixError',NaN,'ComplexTruthError',NaN, ...
    'VoltageRowsError',NaN,'ExactVoltageRankCertificate',false,'SeedRegenerationExact',false,'MatrixKey',0, ...
    'NumericalRank',NaN,'QOrthogonalityError',NaN,'QRangeError',NaN,'RTriangularError',NaN, ...
    'RReciprocalCondition',NaN,'ProjectorIdempotenceError',NaN,'ProjectorAdjointError',NaN, ...
    'IndependentProjectorError',NaN,'ProjectedObservationError',NaN,'OracleAgreementError',NaN, ...
    'MonotonicityMinimum',NaN,'LipschitzMaximumRatio',NaN,'OrthantProjectionError',NaN, ...
    'CircuitMetadataAvailable',false,'CircuitRowsChecked',false,'CircuitRowsError',NaN, ...
    'Pass',false,'ErrorIdentifier','','ErrorMessage','');
end
function need(condition,message)
if ~(isscalar(condition)&&logical(condition)),error('PMUDataAudit:Mismatch','%s',message);end
end
function s=push(s,r)
if isempty(s),s=r;else,s(end+1)=r;end
end
function x=finite_max(values)
values=values(isfinite(values));if isempty(values),x=NaN;else,x=max(values);end
end
function x=finite_min(values)
values=values(isfinite(values));if isempty(values),x=NaN;else,x=min(values);end
end
function restore_session(oldrng,oldpath)
rng(oldrng);path(oldpath);
end
function write_json(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(data,'PrettyPrint',true));
end
