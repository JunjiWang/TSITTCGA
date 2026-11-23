function snap=baseline_snapshot(outfile)
%BASELINE_SNAPSHOT  Record NF/Itr/NormF for a fixed subset of the C group.
%   Used to prove that removing the exact-root shortcut does not change any
%   number: run before and after the edit and compare the two CSV files.
%   The call, the initial points and the parameter sets follow main.m.
here=fileparts(mfilename('fullpath')); addpath(here,'-begin');
names={'TSITTCGA','AT2DFP','ADCGPM','TSIM-DFPM-PRP','HTTCGP'};
solvers={@TSITTCGA,@AT2DFP,@ADCGPM,@TSIM_DFPM_PRP,@HTTCGP};
wjj=struct('Itr_max',2000,'sigma',0.8,'mu',0.6,'kappa',0.5, ...
    'rho_ls',0.8,'a',0.1,'varsigma',2,'rho',1.98,'eta',1e-3,'tbar',0.05, ...
    'c',0.2,'r',0.05,'epsilon_seq',@(j)1/(j+1)^2);
pAT2=struct('Itr_max',2000,'sigma',1e-4,'varsigma',.5,'rho_ls',.4,'gamma',1.99,'tau',.99,'mu',.01,'chibar',.65,'phi',.25);
pADC=struct('Itr_max',2000,'ell',.5,'sigma',.01,'shrink',.55,'rho1',10,'rho2',50,'nu',.8);
pTSIM=struct('Itr_max',2000,'gamma',.37,'sigma',1e-4,'tau',.5,'rho',.5,'eta',.3);
pHTTCGP=struct('Itr_max',2000,'tbar',.3,'mu',.2,'varsigma',1,'rho',.8,'sigma',1e-2,'nu',.8);
params={wjj,pAT2,pADC,pTSIM,pHTTCGP};
rows={};
rng(2016,'twister');
for id=[1 2 5 8]
    for n=[1000 10000]
        PR=struct('name',sprintf('P%d',id),'n',n,'f',@(x,mode) problems(x,mode,id));
        xinit=4*rand(n,1);
        if id==3, xinit=problems(xinit,2,3); end
        for s=1:numel(solvers)
            try
                [~,nf,it,nfres,info]=solvers{s}(PR,xinit,params{s});
                rows(end+1,:)={sprintf('P%d',id),n,names{s},nf,it,nfres,double(info.success)}; %#ok<AGROW>
            catch err
                rows(end+1,:)={sprintf('P%d',id),n,names{s},NaN,NaN,NaN,0}; %#ok<AGROW>
                fprintf('  P%d n=%d %s -> %s\n',id,n,names{s},err.identifier);
            end
        end
    end
end
T=cell2table(rows,'VariableNames',{'Problem','n','Method','NF','Itr','NormF','Success'});
writetable(T,outfile);
fprintf('已写出 %s（%d 行）\n',outfile,height(T));
snap=T;
end
