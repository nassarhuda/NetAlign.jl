# NetworkAlignment.jl
Network Alignment algorithms 
This is an ongoing project.
As a quick sample run:
```
using MatrixNetworks
using Random
n = 50
p = 0.5
A = sparse(erdos_renyi_undirected(n,p));
r = randperm(50);
B = copy(A);
B = A[r,r];

nG = size(A,1)
nH = size(B,1)
w = ones(nG*nH)./(nG*nH)
maxiter = 10
tol = 1e-12
β = 1.0
Xtame = TAME(A,B,w,β,nG,nH,maxiter,tol)

S,w,li,lj = netalign_setup(A,B,L)
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj)
ma,mb = edge_list(bipartite_matching(xbest,li,lj))
```
current supported algorithms:
* netalignmr
* netalignbp
* isornak
* TAME
* NSD
* EigenAlign

Type `?<algorithm_name>` for more info

Up next:
* LowRankEigenAlign