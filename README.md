
## Aircraft Landing Optimization for Enhanced Passenger Comfort

## OVERVIEW 

This project focuses on optimizing the aircraft landing process to enhance passenger comfort by reducing accelerations during landing. Utilizing advanced modeling, simulation, and constrained numerical optimization, the project aims to minimize vibrations and ensure a smooth touchdown. The project combines flight and ground dynamics into a comprehensive simulation model, allowing for detailed analysis and optimization of the landing procedure.

### Features
- **Flight Model Simulation**: Simulates aircraft dynamics during the descent phase, accounting for aerodynamic forces, gravitational forces, and propulsion effects.
- **Ground Model Simulation**: Models the aircraft's behavior post-touchdown, including interactions with the runway and active suspension systems.
- **Hybrid Automaton**: Seamlessly switches between flight and ground dynamics to accurately model the complete landing process.
- **Simulation Results**: Provides insights into the aircraft's trajectory, vertical acceleration, and optimal control inputs through numerical integration techniques.
- **Optimization Framework**: Employs a Finite Horizon Optimal Control Problem (FHOCP) approach to minimize the uncomfortable accelerations passengers experience during landing.


### Technologies
- **Optimization Algorithms**: Implements both unconstrained and constrained optimization methods to find the optimal control strategies.
- **Numerical Methods**: Employs variable and fixed-step numerical integration methods to simulate the system dynamics.

## Notes
The document report.pdf contains a detailed analysis of the project


## CODES

The main folder is divided in three sub-folders:
-Simulation_new
-Full unconstrained final
-Constrained

### Simulation_new 

This sub-folder contains the scripts to perform analysis on how the switching model for the landing procedure behaves.
To use it simply open the script named fly_ground_simulation.m in MATLAB and run it.
### Full unconstrained final 


The scripts in this sub-folder produce the results for the unconstrained version of the optimization problem.
To use it, simply run the script named Full_unc_main.m

### Constrained 

This sub-folder contains the scripts producing the final results deriving from the constrained optimizaton problem.
To run the complete simulation just open the script named final_main.m in MATLAB and run it.
You may also want to run separate simulations for the ground phase and the flight phase; this is possible by running the scripts
new_main_gr.m and new_main_fl.m respectively.

