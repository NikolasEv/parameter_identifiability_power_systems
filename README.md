# Parameter identifiability in power systems

Code for the data-driven (manifold learning) results of

> N. Evangelou, A. M. Stankovic, I. G. Kevrekidis, M. K. Transtrum,
> **Comparing analytic and data-driven approaches to parameter identifiability: A power systems case study.**
> arXiv: [2412.18663](https://arxiv.org/abs/2412.18663)
>
> *The reference will be updated upon acceptance of the paper in a journal.*

The model is the infinite-bus synchronous generator. Trajectories are simulated in MATLAB for 10,000
random perturbations (up to 10%) of the independent parameters. Output-Informed Diffusion Maps,
Double-DMaps Geometric Harmonics and the inverse function theorem are then used in Python to find which
parameters are identifiable.

## Repository layout

```
Codes/
  Matlab/                          data generation (ode45)
    Main_Updated_10_25_2023_10k.m    <- run this one (10,000 trajectories, writes Data/Model_output.mat)
    Main_Updated_10_25_2023.m        original script (50,000 trajectories, different output name)
    SG_model_infinite_bus.m          model equations
    perturb_param.m                  random parameter perturbation
    plot_perturb_param_results.m     plots the simulated trajectories
  4_Power_Systems/                 manifold learning (Jupyter)
    1_Diffusion_Maps_GH_Computations_Deltas.ipynb
    2_Geometric_Harmonics_Deltas.ipynb
    dm/                              diffusion maps / geometric harmonics package used by notebook 2
Data/
  Phis_Parameters.npy              the six non-harmonic Diffusion Maps coordinates (output of notebook 1)
  Model_output.mat                 not included, generated in step 1
Figures/                           figures are written here
```

Keep this layout. The scripts use relative paths (`../../Data`, `../../Figures`) and must be run from their own folders.

## Requirements

- **MATLAB** (developed with R2023b; only base MATLAB is needed) for step 1.
- **Python 3.11** with the packages in `environment.yml`:
  ```bash
  conda env create -f environment.yml
  conda activate param-ident-ps
  ```
  The versions are pinned on purpose. The `dm` package needs `scipy<1.14`.
- **A LaTeX installation** (`latex`, `dvipng`, `ghostscript`), because the figures are drawn with
  `matplotlib`'s `usetex=True`, e.g. TeX Live or MacTeX.

## How to reproduce

### 1. Generate the data (MATLAB)

```matlab
cd Codes/Matlab
Main_Updated_10_25_2023_10k
```

This integrates the model with `ode45` (`RelTol = AbsTol = 1e-7`) up to T = 5 with dt = 0.02 for
10,000 random parameter sets, and saves `Data/Model_output.mat`
(`Behaviors`: 10000 x 251 x 6, `Parameter_values`: 10000 x 12).

### 2. Diffusion Maps (notebook 1)

Open `Codes/4_Power_Systems/1_Diffusion_Maps_GH_Computations_Deltas.ipynb` and run all cells.
This computes the Output-Informed Diffusion Maps and the local linear regression, which finds the six
non-harmonic coordinates. (~20 min on a laptop)

`Data/Phis_Parameters.npy` is already included. To regenerate it from this notebook, uncomment the
`np.save(...)` line in the cell right after `Parameters[:,4] = 0`.

### 3. Geometric Harmonics and inverse function theorem (notebook 2)

Open `Codes/4_Power_Systems/2_Geometric_Harmonics_Deltas.ipynb` and run all cells. It loads
`Data/Model_output.mat` and `Data/Phis_Parameters.npy`. (~20 min on a laptop)

You can also run the notebooks without opening them:

```bash
cd Codes/4_Power_Systems
jupyter nbconvert --to notebook --execute --output-dir /tmp 1_Diffusion_Maps_GH_Computations_Deltas.ipynb
jupyter nbconvert --to notebook --execute --output-dir /tmp 2_Geometric_Harmonics_Deltas.ipynb
```

## Figures and tables

| Paper | File in `Figures/` | Produced by |
|---|---|---|
| Fig. 8  | `Power_System_Local_Linear`             | notebook 1 |
| Fig. 9  | `Power_System_GH_Box_Plots_Parameters`  | notebook 2 |
| Fig. 10 | `Power_Systems_Jacobian_Phis`           | notebook 2 |
| Fig. 11 | `Power_System_GH_Box_Plots_Phis`        | notebook 2 |
| Fig. 12 | `Power_Systems_Jacobian_Parameters`     | notebook 2 |
| Fig. 13 | `Regression_Phis_`                      | notebook 2 |
| Fig. 14 | `Power_Systems_GH_Parity_Plots_Phis`    | notebook 2 |
| Table VI (arXiv v1) | printed output ("Mean Absolute Error Per Parameter") | notebook 2 |

## Note on reproducibility

`perturb_param.m` draws the parameter perturbations with MATLAB's `rand` and no fixed seed, so every
run of step 1 produces a different (statistically equivalent) dataset. Results will then be close to,
but not identical to, those in the paper. With the original dataset, the Python notebooks reproduce the
paper figures exactly. Depending on the BLAS library of your installation, printed values may
differ only at the level of floating-point round-off (~1e-12). If you need that dataset, please contact Nikolas Evangelou through
[https://nikolas-ev.com/](https://nikolas-ev.com/).

## Citation

To be updated upon acceptance of the paper in a journal. For now, please cite the arXiv preprint:

```bibtex
@misc{evangelou2024comparinganalyticdatadrivenapproaches,
      title={Comparing analytic and data-driven approaches to parameter identifiability: A power systems case study}, 
      author={Nikolaos Evangelou and Alexander M. Stankovic and Ioannis G. Kevrekidis and Mark K. Transtrum},
      year={2024},
      eprint={2412.18663},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2412.18663}, 
}
```
