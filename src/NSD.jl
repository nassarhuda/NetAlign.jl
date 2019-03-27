"""

Based on [https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2737&context=cstech]
NSD(A,B,z,w,alpha,iters,m,n)
A,B are the adjacency matrices of the two graphs
z,w are the low rank vectors of the similarity matrix
alpha akin to teleportation parameter
iters is max iters
"""


"""
Network Similarity Decomposition (NSD)
-------
    Based on [https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2737&context=cstech]

    Usage
    ----
    X = NSD(A,B,alpha,maxit,Z,W)
    Input:
    ----
    - `A`: Adjacency matrix of the first graph
    - `B`: Adjacency matrix of the first graph
    - `alpha`: same as alpha in the PageRank equation
    - `maxit`: maximum number of iterations 
    - `Z`: Vector or matrix of the low rank left factor of a similarity matrix L
    - `W`: Vector or matrix of the low rank right factor of a similarity matrix L
      * Z*W' == L

    Methods:
    -------
    X = NSD(A,B,alpha,maxit,z,w) # z and w are vectors
    X = NSD(A,B,alpha,maxit,Z,W) # z and w are matrices
    U,V = NSD_lowrnak(A,B,alpha,maxit,z,w) # z and w are vectors
    U,V = NSD_lowrank(A,B,alpha,maxit,Z,W) # z and w are matrices

    Example:
    -------
    X = NSD(A,B,0.8,10,ones(size(A,1)),ones(size(A,1)))
    ma,mb = edge_list(bipartite_matching(sparse(X))
"""
function NSD(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int},alpha,iters,Zvecs::Array{Float64,2},Wvecs::Array{Float64,2})
  A = _normout_rowstochastic(A)
  B = _normout_rowstochastic(B)
  nB = size(B,1)
  nA = size(A,1)


  # A and B are now row stochastic, so no need for A'x or B'x anywhere
  # operations needed are only Ax or Bx
  GlobalSim = zeros(nA,nB)
  for i = 1:size(Zvecs,2)
    z = Zvecs[:,i]
    w = Wvecs[:,i]
    z./=sum(z)
    w./=sum(w)


    Z = zeros(nA,iters+1)#A
    W = zeros(nB,iters+1)#B
    Sim = zeros(nA,nB)

    W[:,1] = w
    Z[:,1] = z

    for k = 2:iters+1
      W[:,k] = B'*W[:,k-1]
      Z[:,k] = A'*Z[:,k-1]
    end

    for k = 1:iters
      Sim = Sim+alpha^(k-1)*Z[:,k]*W[:,k]'
    end
    Sim = (1-alpha)*Sim
    Sim = Sim + alpha^(iters)*Z[:,iters+1]*W[:,iters+1]'

    GlobalSim += Sim
  end
  return GlobalSim
end

function NSD_lowrank(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int},alpha,iters,Zvecs::Array{Float64,2},Wvecs::Array{Float64,2})
  A = _normout_rowstochastic(A)
  B = _normout_rowstochastic(B)
  nB = size(B,1)
  nA = size(A,1)

  # A and B are now row stochastic, so no need for A'x or B'x anywhere
  # operations needed are only Ax or Bx
  GlobalFactorA = zeros(nA,0)
  GlobalFactorB = zeros(nB,0)
  for i = 1:size(Zvecs,2)
    z = Zvecs[:,i]
    w = Wvecs[:,i]
    z./=sum(z)
    w./=sum(w)


    Z = zeros(nA,iters+1)#A
    W = zeros(nB,iters+1)#B
    Sim = zeros(nA,nB)

    W[:,1] = w
    Z[:,1] = z

    for k = 2:iters+1
      W[:,k] = B'*W[:,k-1]
      Z[:,k] = A'*Z[:,k-1]
    end

    for k = 1:iters
      W[:,k] = sqrt((1-alpha)*alpha^(k-1))*W[:,k]
      Z[:,k] = sqrt((1-alpha)*alpha^(k-1))*Z[:,k]
    end

    W[:,iters+1] = sqrt(alpha^(iters))*W[:,iters+1]
    Z[:,iters+1] = sqrt(alpha^(iters))*Z[:,iters+1]

    GlobalFactorA = hcat(GlobalFactorA,Z)
    GlobalFactorB = hcat(GlobalFactorB,W)
  end
  return GlobalFactorA,GlobalFactorB
end
function NSD(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int},alpha,iters,z::Vector{Float64},w::Vector{Float64})
	A = _normout_rowstochastic(A)
  B = _normout_rowstochastic(B)
  m = size(B,1)
  n = size(A,1)
  z./=sum(z)
  w./=sum(w)
  # A and B are now row stochastic, so no need for A'x or B'x anywhere
  # operations needed are only Ax or Bx

  Z = zeros(n,iters+1)#A
  W = zeros(m,iters+1)#B
  Sim = zeros(m,n)

  W[:,1] = w
  Z[:,1] = z

  for k = 2:iters+1
    W[:,k] = B'*W[:,k-1]
    Z[:,k] = A'*Z[:,k-1]
  end

  # old code (like in the original paper)
  #=
  for k = 1:iters
    Sim = Sim+alpha^(k-1)*W[:,k]*Z[:,k]'
  end
  Sim = (1-alpha)*Sim
  Sim = Sim + alpha^(iters)*W[:,iters+1]*Z[:,iters+1]'
  Sim = Sim' #(Now Sim(a,b) =  similarity between node a in A and node b in B)
  =#
  for k = 1:iters
    Sim = Sim+alpha^(k-1)*Z[:,k]*W[:,k]'
  end
  Sim = (1-alpha)*Sim
  Sim = Sim + alpha^(iters)*Z[:,iters+1]*W[:,iters+1]'
  return Sim
end


function NSD_lowrank(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int},alpha,iters,z,w)
  A = _normout_rowstochastic(A)
  B = _normout_rowstochastic(B)
  m = size(B,1)
  n = size(A,1)
  z./=sum(z)
  w./=sum(w)

  # A and B are now row stochastic, so no need for A'x or B'x anywhere
  # operations needed are only Ax or Bx

  Z = zeros(n,iters+1)#A
  W = zeros(m,iters+1)#B
  Sim = zeros(m,n)

  W[:,1] = w
  Z[:,1] = z

  for k = 2:iters+1
    W[:,k] = B'*W[:,k-1]
    Z[:,k] = A'*Z[:,k-1]
  end

  for k = 1:iters
    W[:,k] = sqrt((1-alpha)*alpha^(k-1))*W[:,k]
    Z[:,k] = sqrt((1-alpha)*alpha^(k-1))*Z[:,k]
  end
  W[:,iters+1] = sqrt(alpha^(iters))*W[:,iters+1]
  Z[:,iters+1] = sqrt(alpha^(iters))*Z[:,iters+1]
  return Z,W
end