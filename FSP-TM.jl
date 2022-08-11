# in this file we have the code for the TM steady state FSP

module FSPTM

    export FSP

    using SparseArrays, LinearAlgebra


    """
    FSP_mat: function that output the FSP matrix for a given state space size N and parameters p
    """
    function FSP_mat(p::Vector,N::Int)
        kon, koff, km, kd = p

        # make sure probability is conserved at i = N
        A00 = sparse([j==i ? -(kd*i + kon) : j==i+1 ? j*kd  : 0 for i = 0:N, j=0:N])
        A01 = sparse([j==i ? koff : 0 for i = 0:N, j=0:N])

        A10 = sparse([j==i ? kon : 0 for i = 0:N, j=0:N])
        A11 = sparse([j==i ? -(kd*i + koff + km) : j==i-1 ? km : j==i+1 ? j*kd  : 0 for i = 0:N, j=0:N])

        A = Matrix(SparseMatrixCSC{Float64,Int64}([A00 A01; A10 A11]))

        # Create the SS FSP matrix with designated state is (1,0)
        As = copy(A)
        As[1,:] = ones(size(As)[1])

        return As
    end

    """
    FSP_mat: solve the FSP problem via the steady state FSP method
    """
    function FSP(p::Vector,N::Int)
        SSmat = FSP_mat(p, N)
        # Create RHS of FSP problem
        b = zeros(size(SSmat)[1]); b[1] = 1;
        # Solve the linear algebra problem.
        sol = SSmat\b
        # get the marginal distribution over M
        M_sol = sol[1:N+1]+sol[N+2:end]
        return M_sol
    end

end
