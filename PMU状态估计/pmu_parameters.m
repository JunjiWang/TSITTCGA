function c=pmu_parameters()
% Frozen application protocol. No tuning is performed by the experiment.
c.systems={'ieee14','ieee118','ieee300'};
c.main_systems={'ieee14','ieee118'};   % main program (IEEE300 runs separately)
c.large_system={'ieee300'};
% Six methods: the five conjugate-gradient projection type solvers (AT2DFP
% included) plus the robust Huber-IRLS reference. WLS and QP are accuracy
% references only and are never ranked by timing.
c.names={'TSITTCGA','AT2DFP','ADCGPM','TSIM-DFPM-PRP','HTTCGP','Huber-IRLS'};
c.solvers={@TSITTCGA_precision,@AT2DFP,@ADCGPM,@TSIM_DFPM_PRP,@HTTCGP};
c.lambda=.005;c.target=1e-6;c.feasibility_target=1e-6;
c.budget=20000;c.repetitions=3;c.backend='normal_eq';c.initial='raw';
% All three systems use the same sparse normal-equation projector.
% Its Cholesky factor is prepared once per instance and shared by all
% methods. Common preparation is excluded from reported solve times.
c.backend_by_system=struct('ieee14','normal_eq','ieee118','normal_eq','ieee300','normal_eq');
c.profile='round2opt_uniform_ne_20261113';
% TSITTCGA: tuned on the PMU problem; the ABC benchmarks keep their own set.
c.params{1}=struct('sigma',.60,'mu',.44,'kappa',.40,'rho',1.99, ...
    'rho_ls',.75,'a',.20,'varsigma',1.33,'eta',1.5e-3,'tbar',.50, ...
    'c',.90,'r',.24,'epsilon_seq',@(j).5,'direction_tolerance',1e-8);
% AT2DFP: fixed signal-recovery configuration of the original paper.
c.params{2}=struct('sigma',.01,'varsigma',1,'rho_ls',.4,'gamma',1.99, ...
    'tau',.99,'mu',.01,'chibar',.65,'phi',.25);
% Peer settings are frozen from the current signal-recovery main.m.
c.params{3}=struct('ell',.5,'sigma',.01,'shrink',.82, ...
    'rho1',10,'rho2',50,'mu',.4,'zeta',50,'xi',.01,'a',2,'gamma',1.9,'eps_F',1e-6);
c.params{4}=struct('gamma',.8,'sigma',.01,'tau',.7,'rho',1.999, ...
    'c',.01,'r',.01,'inertial_fraction',.5,'lambda',.0001,'nu',.01, ...
    'lambda1',.1,'lambda2',.1,'lambda3',.001,'epsilon_seq',@(j).5);
c.params{5}=struct('tbar',.3,'mu',2,'varsigma',1,'rho_ls',.5, ...
    'sigma',.01,'lambda',.001,'nu',.8,'gamma',1.6);
for k=1:5
    c.params{k}.Itr_max=c.budget;c.params{k}.epsilon=c.target;
    c.params{k}.epsilon_C=c.feasibility_target;c.params{k}.quiet=true;
end
end
