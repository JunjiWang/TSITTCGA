function main(options)
%MAIN  Group C: five complete algorithms on the eight monotone benchmarks.
%
%   main                              solve every instance, rebuild tables and figure
%   main(struct('nrep',1))            quick run with one start per instance
%   main(struct('figs_only',true))    redraw the figure only
%   main(struct('tables_only',true))  rebuild the CSV tables only
%
% The work is organised in four stages
%   1. parameters of the five solvers
%   2. test instances and initial points
%   3. solve loop with an independent success check
%   4. persistence, tables and figure

if nargin==0, options=struct; end
o=struct('dims',[10000 30000 50000 80000 100000], ...
    'nrep',5,'Itr_max',2000,'warmup',true,'plots',true, ...
    'figs_only',false,'tables_only',false,'threads',1);
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'Unknown option: %s',key{1});
    o.(key{1})=options.(key{1});
end
validateattributes(o.dims,{'numeric'},{'vector','integer','>=',9,'finite'});
validateattributes(o.nrep,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(o.Itr_max,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(o.threads,{'numeric'},{'scalar','integer','positive','finite'});

% ---- session and folders -----------------------------------------------
here=fileparts(mfilename('fullpath'));
oldpath=path; oldrng=rng; oldfolder=pwd; oldthreads=maxNumCompThreads;
cleanup=onCleanup(@()restore_session(oldpath,oldrng,oldfolder,oldthreads)); %#ok<NASGU>
maxNumCompThreads(o.threads);
cd(here); addpath(here,'-begin');
outdir=fullfile(here,'results'); if ~isfolder(outdir), mkdir(outdir); end
if ~isfolder(fullfile(here,'figures')), mkdir(fullfile(here,'figures')); end
if o.tables_only, make_tables(outdir); return; end
if o.figs_only, plot_results(outdir); return; end

% ---- 1. parameters -----------------------------------------------------
% The proposed method (Algorithm 1 of the manuscript).  This literal
% assignment is parsed verbatim by the Group A/B driver, so it must stay a
% single statement.  The competitor parameters are the ones recommended in
% their own papers (see the header of each solver file).
wjj=struct('Itr_max',o.Itr_max,'sigma',0.8,'mu',0.6,'kappa',0.5, ...
    'a',0.6,'rho_ls',0.8,'varsigma',4,'rho',1.98,'eta',0.01, ...
    'tbar',0.02,'c',0.2,'r',0.15,'epsilon_seq',@(j)1/(j+1)^2,'quiet',true);
pAT2=struct('Itr_max',o.Itr_max,'sigma',1e-4,'varsigma',.5,'rho_ls',.4, ...
    'gamma',1.99,'tau',.99,'mu',.01,'chibar',.65,'phi',.25);
pADC=struct('Itr_max',o.Itr_max,'ell',.5,'sigma',.01,'shrink',.55, ...
    'rho1',10,'rho2',50,'mu',.4,'zeta',50,'xi',.01,'a',2,'gamma',1.9);
pTSIM=struct('Itr_max',o.Itr_max,'gamma',.37,'sigma',1e-4,'tau',.5, ...
    'rho',1.999,'c',.01,'r',.01,'inertial_fraction',.5, ...
    'lambda',.001,'nu',.8,'lambda1',.1,'lambda2',.1,'lambda3',.76);
pHTTCGP=struct('Itr_max',o.Itr_max,'tbar',.3,'mu',.2,'varsigma',1, ...
    'rho_ls',.5,'sigma',.01,'lambda',.001,'nu',.8,'gamma',1.6);
params={wjj,pAT2,pADC,pTSIM,pHTTCGP};
solvers={@TSITTCGA,@AT2DFP,@ADCGPM,@TSIM_DFPM_PRP,@HTTCGP};
names={'TSITTCGA','AT2DFP','ADCGPM','TSIM-DFPM-PRP','HTTCGP'};
nm=numel(solvers);

% ---- 2. instances and initial points -----------------------------------
inst=struct('name',{},'id',{},'n',{},'f',{});
for id=1:8
    for n=o.dims(:)'
        inst(end+1)=struct('name',sprintf('P%d',id),'id',id,'n',n, ...
            'f',@(x,mode)problems(x,mode,id)); %#ok<AGROW>
    end
end
% Every start is 4*rand(n,1).  Problem 3 requires sum(x)<=n, so its start is
% projected onto C; the random stream is unchanged for the other problems.
rng(2016,'twister'); np=numel(inst); X0=cell(np,o.nrep);
for i=1:np
    for r=1:o.nrep
        xinit=4*rand(inst(i).n,1);
        if inst(i).id==3, xinit=problems(xinit,2,3); end
        X0{i,r}=xinit;
    end
end

% ---- 3. solve loop -----------------------------------------------------
source_files=record_sources(solvers,here);
source_snapshot=cellfun(@(f)fileread(fullfile(here,f)),source_files,'UniformOutput',false);
if o.warmup
    for i=1:np
        for j=1:nm
            warm=params{j}; warm.Itr_max=min(3,o.Itr_max);
            solvers{j}(inst(i),X0{i,1},warm);
        end
    end
end
Itr=nan(np,o.nrep,nm); NF=Itr; CPU=Itr; G=Itr; Feasibility=Itr;
RawNF=Itr; RawTime=Itr; Status=cell(size(Itr)); SolutionHash=Status;
Success=false(size(Itr));
timer=tic;
for i=1:np
    for r=1:o.nrep
        for j=mod(r-1+(0:nm-1),nm)+1   % rotate the solver order over starts
            [tc,nf,it,g,d]=solvers{j}(inst(i),X0{i,r},params{j});
            fx=inst(i).f(d.x,1); px=inst(i).f(d.x,2);
            checked=norm(fx); feasible=norm(d.x-px);
            good=d.success && isreal(d.x) && all(isfinite(d.x)) && isreal(fx) ...
                && isfinite(checked) && checked<=1e-6 && isfinite(feasible) && feasible<=1e-6 ...
                && all(isfinite([tc,nf,it,g])) && it<=o.Itr_max;
            assert(~d.success || good,'Solver reported an invalid success.');
            RawNF(i,r,j)=d.actual_NF; RawTime(i,r,j)=d.actual_time;
            Status{i,r,j}=d.status; Success(i,r,j)=good;
            SolutionHash{i,r,j}=sha256(d.x);
            G(i,r,j)=checked; Feasibility(i,r,j)=feasible;
            if good, Itr(i,r,j)=it; NF(i,r,j)=nf; CPU(i,r,j)=tc; end
        end
    end
    fprintf('%2d/%d  %s, n=%d completed\n',i,np,inst(i).name,inst(i).n);
end
wall=toc(timer); matlab_version=version; matlab_threads=maxNumCompThreads;

% ---- 4. persistence ----------------------------------------------------
% Keep the previous results instead of silently overwriting them.
if isfile(fullfile(outdir,'compare.mat'))
    stamp=datestr(now,'yyyymmdd_HHMMSS');
    backup=fullfile(outdir,['archive_' stamp]);
    if ~isfolder(backup), mkdir(backup); end
    copyfile(fullfile(outdir,'*.mat'),backup);
    copyfile(fullfile(outdir,'*.csv'),backup);
    fprintf('Previous results archived in %s\n',backup);
end
save(fullfile(outdir,'compare.mat'),'Itr','NF','CPU','G','Feasibility','Success', ...
    'RawNF','RawTime','Status','SolutionHash','names','params','wjj','o','inst','X0', ...
    'source_files','source_snapshot','wall','matlab_version','matlab_threads','-v7');
make_tables(outdir);
if o.plots, plot_results(outdir); end
fprintf('Saved %d runs. Formal solve time %.2f s.\n',numel(Itr),wall);
end

% =======================================================================
function files=record_sources(solvers,here)
% main.m, problems.m, make_tables.m and every solver file actually used.
files={'main.m','problems.m','make_tables.m'};
for j=1:numel(solvers)
    info=functions(solvers{j});
    assert(strcmpi(fileparts(info.file),here),'Wrong solver source: %s',info.file);
    [~,base,ext]=fileparts(info.file);
    files{end+1}=[base,ext]; %#ok<AGROW>
end
files=unique(files,'stable');
end

function h=sha256(x)
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(double(x(:)),'int8'));
h=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end

function restore_session(savedpath,savedrng,savedfolder,savedthreads)
path(savedpath); rng(savedrng); cd(savedfolder); maxNumCompThreads(savedthreads);
end
