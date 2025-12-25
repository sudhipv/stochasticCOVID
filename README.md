# Stochastic COVID-19 SEIRD Modeling

This repository collects the FreeFEM++, MATLAB, and supporting assets that Sudhi Sharma P V used to study stochastic behavior in a spatially distributed SEIRD (Susceptible–Exposed–Infected–Recovered–Deceased) model for COVID‑19.  
Two complementary uncertainty-quantification strategies live side-by-side:

- **Non-intrusive Monte Carlo** simulations of a square domain with log-normal diffusion fields.
- **Intrusive stochastic spectral finite element (SSFEM)** solvers that propagate random inputs directly in the governing equations for both canonical square problems and a calibrated Southern Ontario geography.

The codebase is intentionally modular so you can re-run the original studies, use individual components (e.g., tensor generators, plotting utilities), or plug in new epidemiological assumptions.

## Highlights
- SEIRD reaction–diffusion solver written in FreeFEM++ with PETSc-backed linear algebra and domain decomposition helper macros.
- Log-normal diffusion sampling (`MonteCarlo/covid_square/mcs_cd.m`) feeding thousands of stochastic runs (`covid_mcs.edp`) whose outputs can be re-aggregated and exported as VTK fields.
- Intrusive SSFEM implementations for one and two random variables, including mesh-specific initializations, higher-order polynomial chaos expansions, and pre-computed Cijk/Tijkl multiplication tensors.
- GIS-informed Southern Ontario case study that blends PHU-level data, calibrated transmission coefficients, and Gaussian reconstructions of the initial condition across 34 health units.
- MATLAB notebooks/scripts that benchmark the intrusive approach against Monte Carlo for SIR/ODE analogs and provide reusable routines for generating chaos tensors.

## Repository Layout
| Path | Purpose |
| --- | --- |
| `MonteCarlo/covid_square/covid_mcs.edp` | Main FreeFEM++ Monte Carlo solver for a square domain, reading `nu.txt` diffusion samples and writing compartment-wise snapshots. |
| `MonteCarlo/covid_square/mcs_cd.m` | MATLAB script that draws Hermite-chaos samples of two log-normal diffusion coefficients and saves them to `nu.txt`. |
| `MonteCarlo/covid_square/plotting/` | MATLAB and FreeFEM++ post-processing utilities (averages, pointwise PDFs, node histories, and VTK exports). |
| `intrusive/covid_1RV/` | Intrusive SSFEM (scalar random variable) for a square domain; contains solver (`covid_ssNL.edp`), template generator (`generate_covid.edp`), weak-form fragments, and 1D chaos tensors (`Tijkl/`, `cijk_O10_D1`). |
| `intrusive/covid_2RV_square/` | Two-RV intrusive solver on a square mesh, plus helper scripts (`generate_2rv.edp`, `uknsum.edp`) and aggregated Monte Carlo comparisons (`compartment*.dat`). |
| `intrusive/covid_2RV_southON/` | Southern Ontario intrusive model with gmsh meshes, PHU data (`PHU*.dat`), Gaussian amplitude tables, and visualization helper `plot_south.m`. |
| `intrusive/MATLAB/SIR_model/` | Deterministic, Monte Carlo, and intrusive SIR ODE solvers used to validate the stochastic numerics (includes `Tijkl_O10_D1`, normalization factors, and comparison figures). |
| `intrusive/MATLAB/TimeDependant/` | Simple ODE benchmark (`ODE_intrusive.m`, `ODE_mcs.m`) plus tensor generators for time-dependent toy problems. |
| `intrusive/MATLAB/NL_spring_matlab/` | Hermite moment/tensor generation scripts (e.g., `Cijk_HPC.m`, `Tijkl_HPC.m`) for the nonlinear spring test case referenced in the thesis. |
| `intrusive/mesh/` & `intrusive/misc/` | Canonical gmsh meshes (various resolutions of squares), pre-generated chaos tensors (`cijk_10_1`), and legacy solver variants. |

## Dependencies
- **FreeFEM++** (with PETSc, MPI, `macro_ddm.idp`, `iovtk`, and `gmsh` modules available in your installation).
- **MPI launcher** such as `mpirun` / `ff-mpirun` for parallel runs.
- **MATLAB** (tested with R2020a) for sample generation, statistical post-processing, and tensor assembly. Most scripts only rely on core MATLAB; `ksdensity` is part of Statistics Toolbox.
- **UQTk** (Uncertainty Quantification Toolkit) if you need to re-create multi-index tables or Hermite quadrature in the tensor scripts; paths in `Cijk_HPC.m` expect an installed copy.
- **gmsh** to regenerate or edit the meshes under `intrusive/mesh` and `intrusive/covid_2RV_southON/mesh`.

You will also need to create writable folders such as `./results/mcs_covid` and `./results/processed` before running the FreeFEM++ post-processing stages, because the scripts assume these directories already exist.

## Monte Carlo Workflow (Non-intrusive)
1. **Generate diffusion samples**  
   Open MATLAB in `MonteCarlo/covid_square` and run:
   ```matlab
   mcs_cd    % edits allowed: number of samples n, distribution params mug/sig
   ```
   This produces the two-column `nu.txt` file (nu\_S and nu\_I). You can edit `mug_0`, `sig_0`, `mug_1`, `sig_1`, and `n` to match a new study.

2. **Launch the FreeFEM++ solver**  
   Ensure you have the `nu.txt` file and (optionally) tweak runtime parameters via `-ARGV`. Example:
   ```bash
   ff-mpirun -n 8 covid_mcs.edp -T 30 -deltaT 0.1 -NLfact 3 -sn 0 -ss 500 \
             -tolP 1e-6 -m 50 -n 50
   ```
   Key flags:
   - `-T` end time and `-deltaT` size, setting `nt`.
   - `-NLfact` scales the initial total population.
   - `-sn` / `-ss` choose which rows of `nu.txt` to sample.
   - `-m`, `-n` control the square mesh resolution if you enable the commented `square(m,n)` block.
   The script writes compartment averages, nodal histories, and nodal fields into `./results/mcs_covid`. Create that directory beforehand.

3. **Process Monte Carlo outputs**  
   MATLAB utilities under `plotting/` provide multiple aggregation modes:
   - `processdata_covid.m`: stacks time-dependent infected averages from batches of 1000 samples, reporting global mean / SD and pointwise PDFs.
   - `process_node_covid.m`: focuses on selected nodes, comparing MCS to intrusive statistics (`./results/processed/intr_nodeI*.dat`) if present.
   - `process_spatial_covid.m`: builds spatial means/variances for the susceptible field, with kernel-density estimates at chosen nodes.
   After computing statistics, export VTK fields via FreeFEM++:
   ```bash
   ff-mpirun -n 1 plotting/plotmcs_covid.edp
   ```
   This reads `./results/processed/mean_nodespaceS.dat` and `sd_nodespaceS.dat` (ensure MATLAB wrote them) and emits VTK files for visualization.

4. **Compare with intrusive results** (optional)  
   `plotting/processpce_covid.m` reconstructs PDFs from intrusive chaos coefficients located in `./output/u*_intr.csv` and compares them to the Monte Carlo histograms for a selected node.

## Intrusive SSFEM Workflow
The intrusive path is split into three practical targets.

### 1. Square domain with one RV (`intrusive/covid_1RV`)
1. **Generate weak-form templates** – run  
   ```bash
   ff-mpirun -n 1 generate_covid.edp -ordi 2 -ordo 3 -muG -12.444 -sigmaG 0.1
   ```  
   This creates `ssinit.edp`, `ssweakcomp.edp`, and companion files that encode the chosen polynomial chaos orders and log-normal moments. The generator reads the `cijk_O10_D1` tensor and produces `Cd.dat`.
2. **Launch the solver** – execute `covid_ssNL.edp` with your target mesh resolution (`-m`, `-n`), time horizon (`-T`, `-deltaT`), and nonlinear scaling (`-NLfact`). The solver pulls together initial conditions, applies the recursive multigrid assembly from `macro_ddm.idp`, and evolves the five compartments.
3. **Clean generated artifacts** using `sh clear.sh` if you want to rebuild with different orders.

### 2. Square domain with two RVs (`intrusive/covid_2RV_square`)
- `generate_2rv.edp` mirrors the single-RV generator but writes chaos definitions that pair two independent log-normal inputs (diffusion in {S,E,R} and in {I}). Run it before the solver whenever `-ordi`/`-ordo` change.
- `covid_2rv.edp` is the main simulator; it shares the same FreeFEM++ backbone but now expands each field into `(ordo+1)` stochastic modes per compartment.  
  Typical launch:
  ```bash
  ff-mpirun -n 8 covid_2rv.edp -m 50 -n 50 -T 0.5 -deltaT 0.1 -NLfact 3 \
            -ordi 2 -ordo 3 -muG1 -12.4 -sigmaG1 0.1 -muG2 -18.4 -sigmaG2 0.1
  ```
- Support files:
  - `uknsum.edp` defines helper macros for the total population chaos expansion.
  - `compartment.dat` / `compartment_sd.dat` capture mean/standard deviation for each compartment over time, enabling later comparisons with Monte Carlo results.
  - `covid_2rv_intrsvforMCS.edp` exports intrusive statistics in the same layout consumed by the Monte Carlo MATLAB scripts.

### 3. Southern Ontario two-RV case (`intrusive/covid_2RV_southON`)
- Uses gmsh meshes under `mesh/` (e.g., `SouthON.msh`, coarse/fine variations).  
- Reads PHU-level coordinates (`mesh/phu_loc.txt`), reported case counts (`mesh/phu_2020_09_01.txt`), and area weights to rebuild Gaussian bumps per health unit. The `B_area.txt` amplitudes were tuned (see `B_valuesoptimized.txt`) so that each Gaussian integrates to the recorded population.
- Transmission coefficients `betaE`/`betaI` and diffusion values are spatially varying logistic combinations fitted for western/central/eastern Ontario.
- To run:
  ```bash
  ff-mpirun -n 16 covid_2rv_south.edp -T 0.4 -deltaT 0.1 -ordi 2 -ordo 3 \
            -muG1 -12.4 -sigmaG1 0.1 -muG2 -18.4 -sigmaG2 0.1
  ```
- The solver reports PHU-level compartment tallies, writes `compartment*.dat`, and can be paired with `plot_south.m` to visualize temporal trends.

## MATLAB Experiments and Utilities
- **`intrusive/MATLAB/SIR_model/`**  
  - `sir_deter.m`: deterministic implicit Euler baseline.  
  - `sir_mcs.m`: Monte Carlo SIR (15k samples) with Picard linearization and PDF estimation.  
  - `sir_intrusive.m`: intrusive SIR using pre-computed `Tijkl_O10_D1`, showcasing Picard convergence and chaos-based variance estimates.  
  - `norm_squared070001.mat` and `Tijkl_O10_D1` carry normalization constants and quadruple moments used both here and in PDE studies.
- **`intrusive/MATLAB/TimeDependant/`**  
  - `ODE_mcs.m` / `ODE_intrusive.m`: simple du/dt + k(θ) u benchmarks for verifying tensor assembly and time-marching.  
  - `Cijk_HPC.m`: generates Hermite chaos tensors by calling UQTk (`generate_quad`, `gen_mi`).  
  - `shade.m` and the `figs/` folder contain plotting helpers used in the thesis.
- **`intrusive/MATLAB/NL_spring_matlab/`**  
  - Houses reusable routines for evaluating third- and fourth-order Hermite integrals (`Cijk_HPC.m`, `Tijkl_HPC.m`, `Sijklm_HPC.m`) and example problems (linear, quadratic, and cubic spring folders).  
  - These tensors feed back into the FreeFEM++ generators whenever you need higher chaos orders.

## Input Files and Supporting Data
- `nu.txt` (repo root under `MonteCarlo/covid_square`) – default 10k log-normal diffusion samples created by `mcs_cd.m`. Regenerate if you adjust `mug`/`sig`.  
- `PHU*.dat` in `intrusive/covid_2RV_southON/` – compartment time series for each Public Health Unit; used for validation/calibration.  
- `Td.dat`, `Cd.dat`, `Tijkl/` directories – multiplication tensors and time-dependent coefficients shared by the generators.  
- `mesh/*.msh` and `.geo` files – gmsh definitions for square and Southern Ontario domains at multiple resolutions.  
- `compartment*.dat`, `MCS*.dat` – sample outputs retained for documentation; feel free to overwrite when running new experiments.

## Tips
- The FreeFEM++ scripts silence most output unless `debug` (Monte Carlo) or `verbosity` flags are set. Toggle them when diagnosing convergence.
- Many solvers read from relative paths (`./results/...`, `./mesh/...`). Run them from their containing directory or use `-cd` to avoid missing-file errors.
- When changing chaos orders, rerun the corresponding `generate_*.edp` file and remove the cached `ss*.edp` pieces with `clear.sh`; otherwise FreeFEM++ will reuse stale definitions.
- Several MATLAB scripts expect Statistics Toolbox functions (`ksdensity`) and write `.dat` files with the `-ascii` flag. Ensure you have write permission in place.

## Questions?
Contact : Sudhi Sharma P V  
Email: sudhisharmapadillath@gmail.com
