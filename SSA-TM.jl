# in this file we have the code for the telegraph model SSA.

module SSATM

    export SSA

    """
    propensity: function to output the propensities given the state of the system.

    args:
    - n: the current state vector.
    - pars: the system parameters.
    """
    function propensity(n::Vector{Int64},pars::Vector{Float64})
        # Reaction rates
        kon, koff, km, kd = pars;
        f_r = zeros(4); # hard code length

        # molecules numbers from the state vector
        ng, nm = n

        # construct the propensities
        f_r[1] = koff*ng;
        f_r[2] = kon*(1-ng);
        f_r[3] = km*ng;
        f_r[4] = kd*nm;

        return f_r::Vector{Float64}
    end

    """
    SSAdt: function to perform the SSA (with the delayed degradation of nRNA) given the necessary parameters.

    args:
    - S_time: the number of individual simulations in the ensemble (set to 1 for single trajectory).
    - pars: the parameters for the ensemble simulations.
    - tol_time: the total simulation time to run for.
    - sp: the storage time period, i.e., if sp = 1.0 then final state vector stored every 1.0s.

    returns:
    - the state vector at the specified times.
    """
    function SSA(S_time::Int, pars::Vector, tol_time::Real, sp::Real)

        sp <= tol_time || error("The storage time period must be less than the total simulation time!")

        # M = Number of reactions, N = Number of reactants
        M = 4::Int; N=2::Int;

        # Define stoichiometry matrix
        S_mat = zeros(N,M);
        S_mat[1,:] = [-1, 1, 0, 0];
        S_mat[2,:] = [0, 0, 1, -1];

        times = convert(Array{Float64,1},LinRange(tol_time,0.0,floor(Int,tol_time/sp)+1));

        # Define reactants trjatory vector
        n = zeros(N,S_time,length(times));

        for sim in 1:S_time
            n_temp = [1,0]; # start gene in active state with zero mRNA
            T = 0;
            sim_times = copy(times);

            # define counter m for updating storage.
            m = 1;
            while T < tol_time
                # Step 1: Calculate propensity
                f_r = propensity(n_temp,pars); # propensity of each reaction.
                lambda = sum(f_r); # total propensity of any reaction.

                # Step 2: Calculate tau and next_r index using random number genrators
                r1 = rand(2,1);
                tau = (1/lambda)*log(1/r1[1]);
                next_r = findfirst(x -> x>=r1[2]*lambda,cumsum(f_r));

                # Step 3: Until the next reaction update the state vector elements for all inbetween times.
                while T+tau >= sim_times[end]
                    n[1:N,sim,m] = n_temp; # m used here.
                    pop!(sim_times);
                    m += 1;
                    if length(sim_times) == 0
                        break
                    end
                end

                # Step 4: update the system time and state vector
                T += tau;
                prod = S_mat[1:N,next_r];
                for i in 1:N
                    n_temp[i] += prod[i]
                end

            end
            if mod(sim,1000) == 0
                println(sim)
            end
        end
        return n
    end

end
