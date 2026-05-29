# Online EM for T2T Association and Spoofing Detection

Minimal MATLAB reproduction package for the EUSIPCO 2026 paper:

> Online EM for T2T Association and Spoofing Detection Using AIS and Radar Data

The entry point is [`main_experiments_AIS_Radar.m`](main_experiments_AIS_Radar.m). It generates the synthetic AIS/radar track pairs used by the paper experiment, runs the online EM estimator, runs the offline baseline, evaluates association and spoofing-detection performance, and exports the paper-style figures.

## Requirements

- MATLAB R2022b or newer is recommended.
- Statistics and Machine Learning Toolbox: `mvnrnd`, `gamrnd`, `randsample`.
- Sensor Fusion and Tracking Toolbox: `assignkbest`.

The script calls [`assert_prerequisites.m`](assert_prerequisites.m) first, so missing MATLAB functions fail early with a short diagnostic.

## Reproduce

From this folder, run:

```matlab
main_experiments_AIS_Radar
```

or from a shell with MATLAB on `PATH`:

```sh
matlab -batch "main_experiments_AIS_Radar"
```

Outputs are written under `results/`:

- `results/outputdata/T2TResults.mat`
- `results/figures/*.pdf`

The random seed is fixed with `rng(1)` in the main script. The default paper-scale configuration is `T = 5000`, `N = 40`, `K = 2`, and `B = 5`.

## Repository Contents

Only the files needed by `main_experiments_AIS_Radar.m` are intended to be tracked:

- Experiment driver: `main_experiments_AIS_Radar.m`
- Online/offline algorithms: `T2TA_onEM.m`, `T2TA_offEM.m`
- Synthetic track generation: `create_track_pairs.m`, `conditional_sampling.m`, `point_mediator.m`
- Filtering and utilities: `ffpfu.m`, `log_sum_exp.m`, `eval_performance.m`
- Plotting: `plot_results_for_paper.m`, `plot_online_offline_results.m`

Older exploratory scripts, real-data files, generated figures, and `.mat` outputs are intentionally ignored by Git so the public repository stays small and reproducible.

## License

See [`LICENSE`](LICENSE).
