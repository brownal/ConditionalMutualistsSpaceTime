using Distributions

function generate_variable_env(lowFit, highFit, mStateProb, sameStateProb, timeScale, maxTime)

	# Map states to the fitness of infected and uninfected individuals in each patch
    stateToFitness = [
		# MM       MP         PM         PP
		lowFit     lowFit     highFit    highFit; # Patch 1, uninfected
		highFit    highFit    lowFit     lowFit; # Patch 1, infected
		lowFit     highFit    lowFit     highFit; # Patch 2, uninfected
		highFit    lowFit     highFit    lowFit] # Patch 2, infected

    # Fraction of the time the two patches are in each state
    mpFreq = pmFreq = (1 - sameStateProb)/2 # Patch 1 in state M & Patch 2 in state P (mpFreq) or the reverse (pmFeq)
    mmFreq = mStateProb - (1 - sameStateProb)/2 # Both patches in State M
    ppFreq= 1 - mmFreq - mpFreq - pmFreq # Both patches in State P

    patchFreqs = [mmFreq, mpFreq, pmFreq, ppFreq]

    # Check that all of the state frequencies are between 0 and 1
    if any(map(x -> (x < 0) || (x > 1), patchFreqs))
        print("Impossible patch state frequencies. sameStateProb is likely too low.")
        return()
    end
    
    # State exit probabilities (i.e. probability the patches leave a given pair of states)
    exitProbs = 1 ./ (timeScale * patchFreqs)

    # Choose an initial state
    currentState = rand(Categorical(patchFreqs))
    currentTime = 1

    # Record the initial time (first column) and state (second column)
    stateTransitions = [currentTime    currentState]

    # Simulate state changes until the time elapsed = maxTime
    while currentTime <= maxTime

        # Find the time at which the patches transition to a new pair of states
        currentTime += rand(NegativeBinomial(1, exitProbs[currentState])) + 1

        # Stop immediately if currentTime is above maxTime
        if currentTime > maxTime
            break
        end

        # Choose the new pair of states (equal probability of being anything but the current one)
        transitionProbs = (1/3) * [(state != currentState) for state=1:4]
        currentState = rand(Categorical(transitionProbs))

        # Record the new state and the time at which it starts
        stateTransitions = vcat(stateTransitions, [currentTime    currentState])

    end

    # Turn the state transitions into a function
    function getState(patch::Int, symb::Int, time::Int)

        # Find the entry in stateTransitions that corresponds to the input time
        stateIndex = findlast(t -> (t <= time), stateTransitions[:, 1])

        # Get the environmental state at that time
        state = stateTransitions[stateIndex, 2]
        
        # Return the fitness for the input patch and infection status
        # Note that the variable symb is 1 when host is uninfected, 2 when infected
        return(stateToFitness[2 * (patch - 1) + symb, state])
    end

    return(stateTransitions, getState)
end
