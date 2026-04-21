We provide code to reproduce the results of the Numerical Experiments in our work Renewable Expectile Regression for Streaming Data: Estimation and Statistical Inference. Two experimental settings are considered:

* Scenario S1: Fix N = 100, 000 and vary the batch size n to evaluate the impact of batch size on performance.

* Scenario S2: Fix nt = 100 and vary the number of batches b to evaluate the number

All scripts should be executed from the project root directory to ensure that relative paths are correctly resolved.

⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻

Description of Main Simulation Code

* functions.R
    Provides core functions for the ReER algorithm and other competitive methods.
  
* sim.R
    Implements the main simulation procedures to reproduce the results.

⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻⸻

Example Usage
```
Target_tau = 0.5
# ==============================================================================
# 模拟 N = 10000，固定场景
run_experiment(
  sim_time = 500, dgpType = 'fixed',
  simType = "sim1",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau
)
# ==============================================================================

# 模拟 nk = 300，stream场景
run_experiment(
  sim_time = 500, dgpType = 'stream',
  simType = "sim1",
  uDist = "norm",
  doScale = FALSE,
  regFormula = 'Y ~ x1 + x2',
  tau = Target_tau
)
```
