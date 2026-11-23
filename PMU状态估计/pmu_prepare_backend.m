function [inst,meta]=pmu_prepare_backend(inst,backend)
%PMU_PREPARE_BACKEND Prepare equivalent orthogonal-projector implementations.
%   dense_q:   P*v = v-Q*(Q'*v), using the frozen orthonormal range basis.
%   sparse_qr: P*v = v-G*(D\v), using one cached sparse QR decomposition D.
%   normal_eq: P*v = v-G*(R\(R'\(G'*v))), where R'*R=G'*G.
%   All variants expose inst.recover(v)=D\v for full-rank least-squares
%   state recovery. The normal_eq variant forms G'*G once per instance.
%
%   This function prepares common model data. Benchmark code must report
%   preparation separately and must not silently charge it to one method only.
if nargin<2,backend='normal_eq';end
backend=char(backend);
assert(any(strcmp(backend,{'dense_q','sparse_qr','normal_eq'})), ...
    'PMU:Backend','Unknown projector backend.');
clock=tic;
G=sparse(inst.G);[m,nstate]=size(G);
assert(m==inst.m && m>=nstate && all(isfinite(nonzeros(G))), ...
    'PMU:BackendDimensions','Invalid measurement matrix.');
factor_clock=tic;
D=decomposition(G,'qr');
factor_seconds=toc(factor_clock);
inst.G=G;
inst.recover=@(v)D\v;
ginfo=whos('G');dinfo=whos('D');qbytes=0;
switch backend
    case 'dense_q'
        assert(isfield(inst,'Q') && isequal(size(inst.Q),[m nstate]), ...
            'PMU:MissingQ','A frozen orthonormal Q is required for dense_q.');
        Q=inst.Q;qinfo=whos('Q');qbytes=qinfo.bytes;
        inst.proj=@(v)v-Q*(Q'*v);
    case 'normal_eq'
        % Equivalent projector through the normal equations of the sparse
        % measurement matrix.  R is nstate-by-nstate and is formed once here.
        R=chol(G'*G);
        inst.proj=@(v)v-G*(R\(R'\(G'*v)));
        inst.Q=[];
    case 'sparse_qr'
        inst.proj=@(v)v-G*(D\v);
        % The projector no longer needs the explicit frozen Q. The input MAT
        % is never modified, and a caller may still retain its own Q reference.
        inst.Q=[];
end
py_clock=tic;inst.Py=inst.proj(inst.y);py_seconds=toc(py_clock);
assert(all(isfinite(inst.Py)),'PMU:BackendNonfinite','Nonfinite projected data.');
meta=struct('Backend',backend,'MeasurementDimension',m,'StateDimension',nstate, ...
    'MeasurementNonzeros',nnz(G),'SparseMeasurementBytes',ginfo.bytes, ...
    'ExplicitQBytes',qbytes,'FactorObjectReportedBytes',dinfo.bytes, ...
    'RecoveryFactorPreparationSeconds',factor_seconds, ...
    'ProjectedDataSeconds',py_seconds,'TotalPreparationSeconds',toc(clock), ...
    'SameExactMathematicalProjector',true, ...
       'NormalEquationsUsed',strcmp(backend,'normal_eq'), ...
    'TimingIncludesPreparation',false);
inst.backend=meta;
end
