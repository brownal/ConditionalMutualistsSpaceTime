# This file is `run_sim_stochastic_environment.jl`.
# Purpose of this is to investigate how stochasticity in the environmental state 
# impacts symbiont containment. This code is run for a burn-in period to allow 
# equilibration, with athe final ~ 32000 time steps recorded. For comparison with
# the deterministic environment case, we set the environment so that the average
# time spent in each pair of states (MM, MP, PM, PP) sums to 32000, the environmental
# cycle length for the deterministic case.
# Infection frequencies and environmental states are recorded for each simulation. 
# Results plotted in ???
#
# This code is meant to run many simulations in parallel. Run it with
# "julia -p # run_sim_stochastic_environment.jl", where # is the
# number of cores to use.
# To run a single simulation, use `myfun` (below) without `pmap`, or
# use the function `simHostCntrl` from sim_host_control.jl

# Parameters
# - Transmission is a host trait
# - No transmission evolution
# - Population initialized at 50% infection in each patch
# - Population size = 200 (100 in each of 2 patches)
# - Simulation run time = 192000 time steps (should see each state roughly 6 times,
#	roughly equivalent to a deterministic time scale of 32000)

# Add the functions for simulating host ecology and generating environments to every core
@everywhere include("sim_host_control.jl")
@everywhere include("generate_random_environment.jl")

# Create lists of the parameter values for each simulation. The idea is to have
# all possible combinations of parameters for the following:
# - flist = fecundity, slist = establishment probability, mlist = mortality;
#   unlike all the other parameters, the values for the symbiont effects only
#   show up with one "var" (variable) and the other two "const" (constant), e.g.
#   (f = "var", s = "const", m = "const") means fecundity varies in time & space
#   as well as with infection status, while establishment and mortality are constant.
# - TSlist = time scale (period) of environmental change
# - dlist = dispersal rate
# - rlist = identifier for each replicate of the simulation
# - filelist = file to save the simulation results in

@everywhere flist    =        [f   for f=("var", "const", "const") for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere slist    =        [s   for s=("const", "var", "const") for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere mlist    =        [m   for m=("const", "const", "var") for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere TSlist   =        [TS  for fsm=1:3                     for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere offlist  = Float64[off for fsm=1:3                     for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere dlist    = Float64[d   for fsm=1:3                     for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere hlist    = Float64[h   for fsm=1:3                     for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere vlist    = Float64[v   for fsm=1:3                     for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere rlist    =        [r   for fsm=1:3                     for d=(0.005, 0.025, 0.05, 0.25, 0.5) for TS=(32000) for off=(0.125, 0.25, 0.375, 0.5, 1.0) for h=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for v=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) for r=1:5]
@everywhere filelist = ["stochastic_environment/time scale = $(TSlist[i]), offset = $(TSlist[i]*offlist[i])/d = $(dlist[i])/f = $(flist[i]), s = $(slist[i]), m = $(mlist[i])/hinit = $(replace("$(hlist[i])", "." => "-")), vinit = $(replace("$(vlist[i])", "." => "-")), sim_$(rlist[i]).csv" for i=1:length(dlist)]


# This function initializes and runs a single simulation
# Input:
# - d: dispersal rate
# - f: fecundity as function of patch, infection status, and time
# - s: establishment probability as function of patch, infection status, and time
# - m: mortality decrease due to infection (M-patches) or lack of infection (P-patches)
# - myfile: file to save results to
# Output:
# - simulation results saved to myfile

@everywhere function myfun(d::Float64, f::String, s::String, m::String, TS::Int, off::Float64, h::Float64, v::Float64, myfile::String)
	# if the directories that contain myfile do not exist, make them & add a file
	# with the parameters
	if !isdir(dirname(myfile))
		mkpath(dirname(myfile))
	end

	# initialize the population so half the hosts in each patch are infected
	initPop = makePop(50*ones(Int, 2, 2), h, v)

	# create functions for the vital rates
	# fec = fecundity, est = establishment, mort = mortality
	if f == "var"
		(fec_env, fec) = generate_random_env(0.5, 1.0, 0.5, 2 * abs(0.5 - off), 32000, 192000);
		
		# save the list of environmental states
		open(myfile[1:(length(myfile) - 4)] * "_env_fec.csv", "w") do f
			write(f, "state_start_time, patch_1, patch_2 \r\n");
			write(f, ["$(fec_env[i, 1]), $(fec_env[i, 2]), $(fec_env[i, 3]) \r\n"
				for i=1:size(fec_env, 1)]);
		end
		
	else
		fec = (patch::Int, symb::Int, t::Int) -> 1.0;
	end

	if s == "var"
		(est_env, est) = generate_random_env(0.5, 1.0, 0.5, 2 * abs(0.5 - off), 32000, 192000);
		
		# save the list of environmental states
		open(myfile[1:(length(myfile) - 4)] * "_env_est.csv", "w") do f
			write(f, "state_start_time, patch_1, patch_2 \r\n");
			write(f, ["$(est_env[i, 1]), $(est_env[i, 2]), $(est_env[i, 3]) \r\n"
				for i=1:size(est_env, 1)]);
		end
		
	else
		est = (patch::Int, symb::Int, t::Int) -> 1.0;
	end

	if m == "var"
		# Note the different order of the first two arguments, because low mortality is better
		(mort_env, mort) = generate_random_env(1.0, 0.5, 0.5, 2 * abs(0.5 - off), 32000, 192000);
		
		# save the list of environmental states
		open(myfile[1:(length(myfile) - 4)] * "_mort_est.csv", "w") do f
			write(f, "state_start_time, patch_1, patch_2 \r\n");
			write(f, ["$(mort_env[i, 1]), $(mort_env[i, 2]), $(mort_env[i, 3]) \r\n"
				for i=1:size(mort_env, 1)]);
		end
		
	else
		mort = (patch::Int, symb::Int, t::Int) -> 1.0;
	end

	# Allow population to reach ecological equilibrium; record final time step
	writeData(simHostCntrl(initPop, d, fec, est, mort,
		0.0, 			# mutation rate = 0.0 (no mutation)
		0.0, 0.0, # standard deviation of horiz. and vert. transmission probability mutations = 0.0
		0.005, 			# spontaneous infection probability = 0.005
		1, 					# number of potentially infectious contacts = 1
		192000, 	# simulation run time = 192000 time steps (six environmental cycles)
		32, 160016), 	# record every 32nd time step, with a burn-in of 160016 time steps (five environmental cycles plus 16 extra time steps for recordings that fall equally at the start and end of each environmental state)
		myfile, 32, 160016);
end

# Map the function to every parameter value (pmap does this in parallel)
pmap((d, f, s, m, TS, off, h, v, myfile) -> myfun(d, f, s, m, TS, off, h, v, myfile),
	dlist, flist, slist, mlist, TSlist, offlist, hlist, vlist, filelist)
