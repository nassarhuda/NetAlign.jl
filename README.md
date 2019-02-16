# NetAlign.jl
Network Alignment algorithms 
This is an ongoing project.
As a quick sample run:
```
S,w,li,lj = netalign_setup(A,B,L)
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj)
ma,mb = edge_list(bipartite_matching(xbest,li,lj))
```
current supported algorithms:
-netalignmr
-netalignbp
-isornak
-TAME
-NSD

Type `?<algorithm_name>` for more info

Up next:
-EigenAlign
-LowRankEigenAlig
