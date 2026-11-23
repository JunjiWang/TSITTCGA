# TTCGPM — Two-Step Inertial Three-Term CG Projection Method

MATLAB code for the paper

> **A Two-Step Inertial Three-Term CG Projection Method for Constrained Equations
> with Application to PMU State Estimation**

The package contains the proposed solver, the four CG-projection methods used for
comparison, a Huber–IRLS baseline, and the scripts that reproduce every table and
figure of the numerical section.

## Requirements

* MATLAB R2025b (any recent release should work; the code uses only base MATLAB)
* [MATPOWER 7.1](https://matpower.org/) — only needed if you want to rebuild the
  PMU input files from the IEEE cases; the shipped `input_ieee14.mat` already
  allows the IEEE 14-bus experiment to run out of the box.

## Layout

| Folder | Contents |
|---|---|
| `单调方程组/` | Group C: eight monotone test problems, five CG projection methods, performance profiles (Fig. 3, Table 4) |
| `A组_线搜索准则/` | Group A: four acceptance criteria compared (Fig. 1, Table 2) |
| `B组_步长策略/` | Group B: warm start and steplength enlargement ablation (Fig. 2, Table 3) |
| `线搜索对比/` | Driver shared by Groups A and B (`main_line_search.m`) |
| `PMU状态估计/` | ℓ1-robust PMU state estimation: efficiency (Table 6), QP verification (Table 7), accuracy (Table 8), Figs. 4–5 |

Root files: `profile_panel.m` (Dolan–Moré panel), `figure_frame.m` (shared
legend), `run_all.m` (run everything), `check_parity.m` (numeric parity check).

## Quick start

```matlab
% Group C, including the performance profiles of Fig. 3
cd 单调方程组
main

% Groups A and B (redraw from stored results only)
cd ../A组_线搜索准则
main(struct('figs_only',true))

% PMU experiments on IEEE 14 (about one minute)
cd ../PMU状态估计
run_pmu('check')
```

`run_all.m` in the repository root runs the three groups in sequence.

## Solvers

| File | Method | Reference |
|---|---|---|
| `TSITTCGA.m` | proposed two-step inertial three-term CG projection | this paper |
| `AT2DFP.m` | accelerated three-term derivative-free projection | Liu (2026) |
| `ADCGPM.m` | derivative-free CG projection with enhanced line search | Jiang, Luo & Jian (2026) |
| `TSIM_DFPM_PRP.m` | two-step inertial modified DFPM with PRP direction | Ma et al. (2020) |
| `HTTCGP.m` | hybrid three-term CG projection | Yin et al. (2021) |
| `TSITTCGA_precision.m` | proposed solver in the PMU (split-CNE) setting | this paper |
| `huber_irls_same_model.m` | Huber–IRLS baseline on the same model | classical |

All solvers use the same stopping rule: `||F(x_k)|| <= 1e-6` at the current
iterate, or the iteration budget is exhausted. No stopping test is applied at
inertial or trial points.

## Data

The network data and operating points are the IEEE 14-, 118- and 300-bus test
cases distributed with MATPOWER 7.1, originally from the
[Power Systems Test Case Archive](https://labs.ece.uw.edu/pstca/) of the
University of Washington. The PMU measurements themselves are **synthetic**:
voltage and current phasors are computed from the operating points, Gaussian
noise with standard deviation `1e-3` p.u. per real component is added, and a
prescribed number of components is contaminated by gross errors. Every instance
is generated with a fixed seed, so the shipped scripts reproduce the published
numbers exactly.

`input_ieee14.mat` (2 MB) is included. The IEEE 118- and 300-bus measurement
files are 31 MB and 149 MB and are therefore **not** shipped; rebuild them with

```matlab
build_ieee118_input('/path/to/matpower7.1');
build_ieee300_input('/path/to/matpower7.1');
```

placed in `PMU状态估计/`.

## Reproducing the tables

```matlab
cd 单调方程组;      main            % Table 4 and Fig. 3
cd ../线搜索对比;    main_line_search % Tables 2-3 and Figs. 1-2
cd ../PMU状态估计;   main            % Tables 6-8 and Figs. 4-5 (several hours)
```

Results are written to the `results/` subfolder of each group and are not
tracked by git.

## Notes

* Timings are wall-clock times measured with `tic`/`toc` under a single
  computational thread, averaged over repeated runs.
* The projector `Pi_p = I - A_p A_p^+` is applied through one sparse Cholesky
  factorisation of `A_p' A_p`, formed once per instance and timed separately.
* Every solver receives identical initial points and iteration budgets.


## Runtime and safety

| Entry point | Typical runtime | What it writes |
|---|---|---|
| `单调方程组/main` | 10-20 min | `results/`, `figures/` |
| `A组_线搜索准则/main`, `B组_步长策略/main` | 1-2 min (figures only) | `figures/` |
| `PMU状态估计/run_pmu('check')` | about 6 min | `results/` |
| `PMU状态估计/main` | several hours | `results/` |

Nothing is written outside the folder of the corresponding group. If
`results/compare.mat` already exists, the previous results are first copied to
`results/archive_<timestamp>/` and only then replaced, so a rerun never
silently destroys earlier output.

## License

MIT — see `LICENSE`.

## Contact

Junji Wang (junjiwang@example.com) — please replace with the corresponding
author's address before publishing.
