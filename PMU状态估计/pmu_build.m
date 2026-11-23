function inst=pmu_build(c,V,bus_sel,sigma,seed,nbad,cur_sel)
%PMU_BUILD Exact rectangular PMU measurements and a fixed QR projector.
% Voltage phasors are measured at bus_sel; cur_sel indexes terminal currents
% as 2*k-1 (FROM), 2*k (TO) for branch k. A current is included only if its
% terminal bus is selected. sigma is per real/imaginary p.u. component.
% Selected measurements must jointly observe the full voltage state.
% Reduced model: min_e .5*||P(y-e)||^2+lambda*||e||_1,
% P=I-QQ', where Q is an orthonormal basis of range(G).
n=c.n;V=V(:);nl=numel(c.branch.f);
if nargin<7 || isempty(cur_sel),cur_sel=1:2*nl;end
assert(numel(V)==n && all(isfinite(V)) && all(abs(V)>0), ...
    'PMU:State','Invalid true complex voltage state.');
assert(isreal(sigma) && isscalar(sigma) && isfinite(sigma) && sigma>0, ...
    'PMU:Noise','sigma must be positive and finite.');
assert(isscalar(seed) && isfinite(seed) && seed>=0 && seed==fix(seed) ...
    && seed<=2^32-1,'PMU:Seed','Invalid random seed.');
assert(isreal(bus_sel) && all(isfinite(bus_sel(:))) ...
    && all(bus_sel(:)>=1 & bus_sel(:)<=n & bus_sel(:)==fix(bus_sel(:))) ...
    && numel(unique(bus_sel))==numel(bus_sel),'PMU:BusSelection','Invalid bus selection.');
assert(isreal(cur_sel) && all(isfinite(cur_sel(:))) ...
    && all(cur_sel(:)>=1 & cur_sel(:)<=2*nl & cur_sel(:)==fix(cur_sel(:))) ...
    && numel(unique(cur_sel))==numel(cur_sel), ...
    'PMU:CurrentSelection','Invalid terminal-current selection.');
sel=ismember(1:n,bus_sel);H=zeros(0,n);kind={};
for i=find(sel)
    a=zeros(1,n);a(i)=1;H=[H;a];kind{end+1}=sprintf('V%d',i); %#ok<AGROW>
end
for k=1:nl
    f=c.branch.f(k);tt=c.branch.t(k);
    if sel(f) && ismember(2*k-1,cur_sel)
        a=zeros(1,n);a(f)=c.branch.Yff(k);a(tt)=c.branch.Yft(k);
        H=[H;a];kind{end+1}=sprintf('I%d-%d',f,tt); %#ok<AGROW>
    end
    if sel(tt) && ismember(2*k,cur_sel)
        a=zeros(1,n);a(f)=c.branch.Ytf(k);a(tt)=c.branch.Ytt(k);
        H=[H;a];kind{end+1}=sprintf('I%d-%d',tt,f); %#ok<AGROW>
    end
end
nq=size(H,1);m=2*nq;
assert(isscalar(nbad) && isfinite(nbad) && nbad>=0 && nbad==fix(nbad) && nbad<=m, ...
    'PMU:BadDataCount','Invalid number of corrupted real measurement components.');
G=[real(H),-imag(H);imag(H),real(H)];xi=[real(V);imag(V)];
q=H*V;ytrue=[real(q);imag(q)];
rankG=rank(G);assert(rankG==2*n,'PMU:Unobservable','PMU model is not full column rank.');
[Q,~]=qr(G,0);orth_error=norm(Q'*Q-eye(2*n),'fro');
range_error=norm(G-Q*(Q'*G),'fro')/max(1,norm(G,'fro'));
truth_error=norm(G*xi-ytrue,inf);
assert(orth_error<1e-10 && range_error<1e-10 ...
    && truth_error<1e-10*max(1,norm(ytrue,inf)), ...
    'PMU:Projector','Measurement/QR projector consistency check failed.');
% Preserve the previous random draw order and corruption definition.
rng(seed,'twister');bad=zeros(m,1);idx=[];
if nbad>0
    idx=randperm(m,nbad);
    bad(idx)=(20+80*rand(nbad,1)).*sign(randn(nbad,1))*sigma;
end
noise=sigma*randn(m,1);y=ytrue+noise+bad;
inst=struct('G',G,'y',y,'ytrue',ytrue,'true_xi',xi,'bad',bad, ...
    'nq',nq,'m',m,'n',n,'H',H,'kind',{kind},'sigma',sigma, ...
    'seed',seed,'noise',noise,'bad_indices',idx(:),'nbad',nbad, ...
    'bus_selection',bus_sel(:),'current_selection',cur_sel(:), ...
    'rank',rankG,'cond',cond(G),'Q',Q, ...
    'qr_orthogonality_error',orth_error,'qr_range_error',range_error, ...
    'measurement_identity_error',truth_error);
inst.proj=@(z)z-Q*(Q'*z);
end
