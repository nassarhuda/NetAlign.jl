"""
EigenAlign: Spectral Alignment of Networks
-------
    Based on [https://arxiv.org/pdf/1602.04181v1.pdf]

    Usage
    ----
    ma,mb = EigenAlign(A,B)

    Input:
    ----
    - `A`: Adjacency matrix of the first graph
    - `B`: Adjacency matrix of the first graph
    - `s1`: weight of overlaps
    - `s2`: weight of non-informatives
    - `s3`: weight of conflicts

    Methods:
    -------
    ma,mb = EigenAlign(A,B)
    ma,mb = EigenAlign(A,B,s1,s2,s3)
    ma,mb = EigenAlign(A,B,residerror,maxiter)
    ma,mb = EigenAlign(A,B,s1,s2,s3,residerror,maxiter)

    Example:
    -------
    n = 50; p = 0.2;
    A = sparse(erdos_renyi_undirected(n,p));
    r = randperm(n)
    B = A[r,r]
    ma,mb = EigenAlign(A,B)
    normalized_overlap(A,B,ma,mb)
"""
function EigenAlign(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int};residerror::Float64=1e-12,maxiter::Int=100)
  s1,s2,s3 = find_parameters(A,B)
  return EigenAlign(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int},s1,s2,s3;residerror=residerror,maxiter=maxiter)
end
function EigenAlign(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int},s1::Int,s2::Int,s3::Int;residerror::Float64=1e-12,maxiter::Int=100)
  gam1 = s1+s2-2s3
  gam2 = s3-s2
  gam3 = s2

  nA = size(A,1)
  nB = size(B,1)

  Eb = ones(Int,nB,nB)
  Ea = ones(Int,nA,nA)

  X = 1/(nA*nB)*ones(nB,nA)

  resid = 1
  curiter = 1
  while resid>residerror || curiter<=maxiter
    Ax = gam1*B*X*A' + gam2*Eb*X*A' + gam2*B*X*Ea + gam3*Eb*X*Ea # this is equivalent to Mx from original paper
    y = Ax[:]
    x = X[:]
    lam = (x'*y)./(x'*x)
    resid = norm(y - lam*x)
    X = Ax./norm(y,1)
    curiter+=1
  end
  
  Xmat = reshape(X,nB,nA);

  ej,ei = edge_list(bipartite_matching(sparse(Xmat)))

  if curiter == maxiter
    warn("maximum number of iterations reached with residual = $resid ")
  end

  return ei,ej
end
function find_parameters(A::SparseMatrixCSC{Int,Int},B::SparseMatrixCSC{Int,Int})
  nB = size(B,1)
  nA = size(A,1)
  nmatches = sum(A)*sum(B)
  nmismatches = sum(A)*(nB^2 - sum(B)) + sum(B)*(nA^2 - sum(A))
  mygamma = nmatches/nmismatches
  myalpha = (1/mygamma) + 1
  myeps = 0.001
  s1 = myalpha + myeps
  s2 = 1+myeps
  s3 = myeps
  return s1,s2,s3
end