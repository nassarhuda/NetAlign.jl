# evaluations script
"""
measure of the network alignment result

Example:
-------
S,w,li,lj = netalign_setup(A,B,L)
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj)
ma,mb = edge_list(bipartite_matching(xbest,li,lj))
r = normalized_overlap(A,B,ma,mb)
"""
function normalized_overlap(A,B,ma,mb)
	C = A[ma,ma].*B[mb,mb]
	nnz(C)/max(nnz(A),nnz(B))
end

"""
measure of the network alignment result
Precision and Recall given a ground truth alignmnent

Example:
-------
S,w,li,lj = netalign_setup(A,B,L)
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj)
ma,mb = edge_list(bipartite_matching(xbest,li,lj))
precision,recall,correcmatches = ground_truth_measure(ma,mb,maref,mbref)
"""
function ground_truth_measure(ma,mb,maref,mbref)
	Ma = hcat(ma,mb)
	Mref = hcat(maref,mbref)
	ids = findin_rows(Mans,Mref)
	P = length(ids)/length(ma)
	R = length(ids)/length(maref)
	return P,R,ids
end