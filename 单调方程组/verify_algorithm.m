function report=verify_algorithm(outdir)
%VERIFY_ALGORITHM  Check that the code follows Algorithm 1 of the manuscript.
%
% The solver runs with para.debug_dir=true, which records (x_k, varpi_k,
% F(varpi_k), d_k) at every iteration.  The script then recomputes varpi_k
% and d_k from the formulas of the manuscript (no intermediate quantity of
% the solver is reused) and compares the two.  This is a code-vs-formula
% check only; the descent and trust-region properties are not tested here.
%
% Usage:  report = verify_algorithm;            % writes results/algorithm_verification.csv

solver=functions(@TSITTCGA);
here=fileparts(solver.file);
problem_source=functions(@problems);
assert(strcmpi(fileparts(problem_source.file),here), ...
    'verify_algorithm:Source','The solver and problem definitions must come from the same folder.');
if nargin<1 || isempty(outdir), outdir=fullfile(here,'results'); end
if ~isfolder(outdir), mkdir(outdir); end

main_source=fileread(fullfile(here,'main.m'));
parameter_assignment=regexp(main_source,'(?m)^\s*wjj\s*=\s*struct[\s\S]*?;','match','once');
assert(~isempty(parameter_assignment),'verify_algorithm:Parameters', ...
    'Cannot locate the wjj parameter assignment in main.m.');
o=struct('Itr_max',600); %#ok<NASGU>
eval(parameter_assignment);
para=wjj; para.Itr_max=600; para.debug_dir=true; para.quiet=true;
eps_F=1e-6; eps_C=1e-6;

saved_rng=rng; rng_cleanup=onCleanup(@()rng(saved_rng)); %#ok<NASGU>
rng(314159,'twister');
rows=cell(0,6); nfail=0;
for id=1:8
    for n=[100 1000]
        if id==3, x0=2*rand(n,1)-1; else, x0=4*rand(n,1); end
        PR=struct('name',sprintf('P%d',id),'n',n,'f',@(x,mode)problems(x,mode,id));
        [~,~,~,~,info]=TSITTCGA(PR,x0,para);
        h=info.dir_history;
        assert(~isempty(h),'verify_algorithm:EmptyHistory','No direction history for P%d.',id);
        assert(all(arrayfun(@(q)q.k,h)==(1:numel(h))),'verify_algorithm:IndexGap', ...
            'Direction history is not consecutive.');
        X=[h.x]; W=[h.w]; Fv=[h.F]; D=[h.d];
        maxdir=0; maxw=0; nfallback=0;
        for i=1:numel(h)
            j=h(i).j;
            % ---- varpi_k ---------------------------------------------------
            if j==0
                wref=X(:,i);
            else
                ek=para.epsilon_seq(j);
                dx=norm(X(:,i)-W(:,i-1));
                if dx>0, ck=min(para.c,ek/(j^2*dx)); else, ck=para.c; end
                if j==1
                    wref=X(:,i)+ck*(X(:,i)-W(:,i-1));
                else
                    dr=norm(X(:,i-1)-W(:,i-2));
                    if dr>0, rk=min(para.r,ek/(j^2*dr)); else, rk=para.r; end
                    wref=X(:,i)+ck*(X(:,i)-W(:,i-1))+rk*(X(:,i-1)-W(:,i-2));
                end
            end
            Fref=PR.f(wref,1);
            if norm(Fref)<=eps_F && norm(wref-X(:,i))>eps_C
                wref=X(:,i);
                nfallback=nfallback+1;
            end
            maxw=max(maxw,norm(W(:,i)-wref)/max(1,norm(wref)));
            % ---- d_k -------------------------------------------------------
            if j==0
                dref=-Fv(:,i);
            else
                s=W(:,i)-W(:,i-1);
                y=Fv(:,i)-Fv(:,i-1);
                s2=s'*s; y2=y'*y;
                if s2==0 || y2==0
                    dref=-Fv(:,i);
                else
                    Ahat=-(3/5)*(y'*s);
                    p=max(0,para.eta*norm(Fv(:,i))-(y'*s+Ahat)/s2);
                    tau=min(para.tbar,max(0,1-(y'*s)/y2));
                    ytil=y+(para.eta+Ahat/s2+p)*s;
                    mu=max([norm(Fv(:,i-1))^2,para.varsigma*norm(D(:,i-1))*norm(ytil), ...
                        D(:,i-1)'*ytil]);
                    beta=(Fv(:,i)'*ytil)/mu-(ytil'*ytil)*(Fv(:,i)'*D(:,i-1))/mu^2;
                    theta=tau*(Fv(:,i)'*D(:,i-1))/mu;
                    dref=-Fv(:,i)+beta*D(:,i-1)+theta*ytil;
                end
            end
            maxdir=max(maxdir,norm(D(:,i)-dref)/max(1,norm(dref)));
        end
        ok=maxdir<1e-10 && maxw<1e-10;
        if ~ok, nfail=nfail+1; end
        fprintf(['P%d n=%6d itr=%3d  max|d-dref|/|dref|=%.2e  max|w-wref|/|wref|=%.2e  ' ...
            'fallbacks=%d  %s\n'],id,n,numel(h),maxdir,maxw,nfallback,string(ok));
        rows(end+1,:)={sprintf('P%d',id),n,numel(h),maxdir,maxw,nfallback}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'Problem','n','Iterations', ...
    'MaxDirectionError','MaxInertialPointError','FallbackCount'});
writetable(T,fullfile(outdir,'algorithm_verification.csv'),'Encoding','UTF-8');
report=struct('cases',size(T,1),'failures',nfail, ...
    'max_direction_error',max([T.MaxDirectionError]), ...
    'max_inertial_error',max([T.MaxInertialPointError]), ...
    'fallback_count',sum(T.FallbackCount));
fprintf(['\nIMPLEMENTATION_CHECK cases=%d failures=%d ' ...
    'maxDirErr=%.2e maxInertialErr=%.2e\n'],report.cases,nfail, ...
    report.max_direction_error,report.max_inertial_error);
if nfail==0
    fprintf('IMPLEMENTATION_CHECK_PASSED\n');
else
    fprintf('IMPLEMENTATION_CHECK_FAILED\n');
end
end



