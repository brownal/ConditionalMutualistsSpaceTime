Code for "Evolution of symbiont transmission in conditional mutualisms in spatially and temporally variable environments"

Julia code:
  * sim_host_control.jl: contains functions for simulating transmission evolution in a population of hosts in an environment that changes in space and time.
  * run_sim_ecological_dynamics.jl: simulates ecological dynamics of infection in a population of hosts with fixed transmission probabilities.
  * run_sim_transmission_evolution.jl: simulates transmission evolution when the environment changes over intermediate timescales (160 and 80 generations).
  * run_sim_more_state_M.jl: simulates transmission evolution when patches spend more time in State M (where the symbiont is mutualistic) than State P.
  * run_sim_more_state_P.jl: simulates transmission evolution when patches spend more time in State P (where the symbiont is parasitic) than State M.
  * run_sim_fast_env_change.jl: simulates transmission evolution when the environment changes over very short timescales (4 generations).
  * run_sim_slow_env_change.jl: simulates transmission evolution when the environment changes over very long timescales (40000 generations).

Mathematica code:
  * numerically_finding_containment.wls: finds the average symbiont containment for populations of hosts with fixed transmission probabilities, using numerical solutions of an analytical model.