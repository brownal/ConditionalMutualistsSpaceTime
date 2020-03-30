# This file is `run_sim_fast_env_change.jl`.
# Purpose of this is to simulate transmission evolution when the environment changes
# rapidly. Results plotted in figure S5.
#
# This code is meant to run many simulations in parallel. Run it with
# "julia -p # run_sim_fast_env_change.jl", where # is the
# number of cores to use.
# To run a single simulation, use `myfun` (below) without `pmap`, or
# use the function `simHostCntrl` from sim_host_control.jl

# Parameters
# - Transmission is a host trait
# - Population initialized at 50% infection in each patch
# - Population size = 200 (100 in each of 2 patches)
# - Simulation run time = 2*10^6 time steps

# Add the functions for simulating host evolution to every core
@everywhere include("sim_host_control.jl")

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

@everywhere flist = [f for f=("var", "const") for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere slist = [s for s=("const", "const") for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere mlist = [m for m=("const", "var") for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere TSlist = [TS for fsm=1:3 for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere offlist = Float64[off for fsm=1:3 for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere dlist = Float64[d for fsm=1:3 for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere hlist = Float64[h for fsm=1:3 for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere vlist = Float64[v for fsm=1:3 for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere rlist = [r for fsm=1:3 for d=(0.005, 0.05, 0.5) for TS=(800) for off=(0.25, 0.5, 1.0) for h=(0.0, 0.5, 1.0) for v=(0.0, 0.5, 1.0) for r=1:3]
@everywhere filelist = ["fast_env_change/time scale = $(TSlist[i]), offset = $(TSlist[i]*offlist[i])/d = $(dlist[i])/f = $(flist[i]), s = $(slist[i]), m = $(mlist[i])/hinit = $(replace("$(hlist[i])", "." => "-")), vinit = $(replace("$(vlist[i])", "." => "-")), sim_$(rlist[i]).csv" for i=1:length(dlist)]


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
		fec = (patch::Int, symb::Int, t::Int) -> 0.5*ceil((mod(t + (TS*off)*(patch-1) + (TS/2)*(symb-1), TS)/(TS - 1)) - (1/2)) + 0.5;
	else
		fec = (patch::Int, symb::Int, t::Int) -> 1.0;
	end

	if s == "var"
		est = (patch::Int, symb::Int, t::Int) -> 0.5*ceil((mod(t + (TS*off)*(patch-1) + (TS/2)*(symb-1), TS)/(TS - 1)) - (1/2)) + 0.5;
	else
		est = (patch::Int, symb::Int, t::Int) -> 1.0;
	end

	if m == "var"
		mort = (patch::Int, symb::Int, t::Int) -> 0.5*(1 - ceil((mod(t + (TS*off)*(patch-1) + (TS/2)*(symb-1), TS)/(TS - 1)) - (1/2))) + 0.5;
	else
		mort = (patch::Int, symb::Int, t::Int) -> 1.0;
	end

	# Allow population to reach ecological equilibrium (no evolution)
	writeData(simHostCntrl(initPop, d, fec, est, mort,
		0.02, 			# mutation rate = 0.02
		0.05, 0.05, # standard deviation of horiz. and vert. transmission probability mutations = 0.01
		0.005, 			# spontaneous infection probability = 0.002
		1, 					# number of potentially infectious contacts = 1
		16000000, 	# simulation run time = 1.6*10^7 time steps
		20000 + Int(TS/8), Int(TS/16)), 	# record every 20000 + TS/8 time step
		myfile, 20000 + Int(TS/8), Int(TS/16));
end

# Map the function to every parameter value (pmap does this in parallel)
pmap((d, f, s, m, TS, off, h, v, myfile) -> myfun(d, f, s, m, TS, off, h, v, myfile),
	dlist, flist, slist, mlist, TSlist, offlist, hlist, vlist, filelist)
