# bayesFlow

Codes in this repository generate probabilistic forecasts of international migration flows between the 200 most populous countries. 
Computationally intense components of the implementation are contained in the `cluster` directory. 
Codes used to evaluate and summarize the forecast results are in the `local` directory. 

The `cluster` codes are designed to run as batch processes handled by the Slurm management and job scheduling system. 
We exclude `cluster/{forecast, heldout, lib, log}` from version control to avoid saving multiple gigabytes of simulation results. 
