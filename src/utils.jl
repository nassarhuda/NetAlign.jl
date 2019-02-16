"""
creates a column stochastic matrix
"""
function _normout_colstochastic(P::SparseArrays.SparseMatrixCSC{T,Int64}) where T
	n = size(P,1)
	colsums = sum(P,dims=2)
	pi,pj,pv = findnz(P)
	Q = SparseMatrixCSC(P.m,P.n,P.colptr,P.rowval,pv./colsums[pj])
end
"""
creates a row stochastic matrix
"""
function _normout_rowstochastic(P::SparseArrays.SparseMatrixCSC{T,Int64}) where T
	n = size(P,1)
	colsums = sum(P,dims=2)
	pi,pj,pv = findnz(P)
	Q = SparseMatrixCSC(P.m,P.n,P.colptr,P.rowval,pv./colsums[pi])
end

"""
findin_rows(Mans,Mref) returns an indicator vector of the rows in Mans that are ground truth matches
```
example:
julia> a = [1 2; 
            2 3; 
            3 1];
julia> b = [1 2; 
            3 3; 
            2 1];
julia> findin_rows(a,b) # [1]
```
"""
function findin_rows(Mans,Mref)
  a = [Mans[i,:] for i in 1:size(Mans,1)]
  b = [Mref[i,:] for i in 1:size(Mref,1)]
  ind = Vector{Int}(undef,0)
  bset = Set(b)
  @inbounds for (i,ai) in enumerate(a)
    ai in bset && push!(ind, i)
  end
  ind
end