"""
Triangular Alignment (TAME): A Tensor-based Approach for Higher-order Network Alignment
-------
    Based on [https://arxiv.org/pdf/1510.06482.pdf]

    Usage
    ----
    X = TAME(A,B,w,β,nA,nB,maxiter,tol)

    Input:
    ----
    - `A`: Adjacency matrix of the first graph
    - `B`: Adjacency matrix of the first graph
    - `w`: vectorized version of the sequence similarities between nodes from A and B
    - `beta`: shift parameter used in the paper
    - `nA`: number of nodes in graph A
    - `nB`: number of nodes in graph B
    - `maxit`: maximum number of iterations 
    - `tol`: tolerance accuracy of the eigenvalue obtained

    Methods:
    -------
    X = TAME(G,H,w,β,nG,nH,maxiter,tol)
    X = cTAME(G,H,w,β,nG,nH,maxiter,tol)
    X = TAME(G,H)
    X = cTAME(G,H)

    Example:
    -------
    X = TAME(G,H,w,β,nG,nH,maxiter,tol)
    ma,mb = edge_list(bipartite_matching(sparse(X))
"""
function TAME(G::SparseMatrixCSC{Int64,Int64},H::SparseMatrixCSC{Int64,Int64})
	nG = size(G,1)
	nH = size(H,1)
	w = ones(nG*nH)./(nG*nH)
	maxiter = 1000
	tol = 1e-12
	β = 1
	return TAME(G,H,w,β,nG,nH,maxiter,tol)
end
function cTAME(G::SparseMatrixCSC{Int64,Int64},H::SparseMatrixCSC{Int64,Int64})
	nG = size(G,1)
	nH = size(H,1)
	w = ones(nG*nH)./(nG*nH)
	maxiter = 1000
	tol = 1e-12
	β = 1
	return cTAME(G,H,w,β,nG,nH,maxiter,tol)
end
function TAME(G::SparseMatrixCSC{Int64,Int64},H::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},β::Float64,nG::Int,nH::Int,maxiter::Int,tol::Float64)
	k = 0
	w = w./norm(w,1)
	xcur = w
	Xbest = reshape(xcur,nG,nH)
	bestscore = 0
	itercount = 1
	oldlam = 0
	while itercount <= maxiter
		xnew = impTTV(G,H,xcur,nG,nH)
		lam = dot(xcur,xnew)
		xnew .+= β*xcur
		xnew ./= norm(xnew,1)
		Xmat = reshape(xnew,nG,nH)
		curscore = score_fn(Xmat)
		if curscore > bestscore
			bestscore = curscore
			Xbest = Xmat
		end
		xcur = xnew
		itercount += 1
		if itercount == 1
			oldlam = lam
		else
			if lam-oldlam < tol
				break
			end
			oldlam = lam
		end
	end
	return Xbest
end

function cTAME(G::SparseMatrixCSC{Int64,Int64},H::SparseMatrixCSC{Int64,Int64},w::Vector{Float64},β::Float64,nG::Int,nH::Int,maxiter::Int,tol::Float64)
	k = 0
	w = w./norm(w,1)
	W = reshape(w,nG,nH)
	xcur = w
	Xbest = reshape(xcur,nG,nH)
	bestscore = 0
	itercount = 1
	oldlam = 0
	while itercount <= maxiter
		xnew = c_impTTV(G,H,x,nG,nH,W)
		lam = dot(xcur,xnew)
		xnew .+= β*xcur
		xnew ./= norm(xnew,1)
		Xmat = reshape(xnew,nG,nH)
		curscore = score_fn(Xmat)
		if curscore > bestscore
			bestscore = curscore
			Xbest = Xmat
		end
		xcur = xnew
		itercount += 1
		if itercount == 1
			oldlam = lam
		else
			if lam-oldlam < tol
				break
			end
			oldlam = lam
		end
	end
	return Xbest
end

function c_impTTV(G::SparseMatrixCSC{Int64,Int64},H::SparseMatrixCSC{Int64,Int64},x::Vector{Float64},nG::Int,nH::Int,W::Array{Float64,2})
	X = reshape(x,nG,nH)
	Y = similar(X)
	Y .= 0
	for g = 1:nG
		for h = 1:nH
			g_triangles = triangles(G,g) # no symmetry
			h_triangles = triangles(H,h) # no symmetry
			for g_tri in g_triangles
				j,k = g_tri.w,g_tri.x
				for h_tri in h_triangles
					jp,kp = h_tri.w,h_tri.x
					if W[g,h] !=0
						Y[g,h] += X[j,jp]*X[k,kp]+X[j,kp]*X[k,jp]
					else
						Y[g,h] = 0
					end
				end
			end
			Y[g,h] = Y[g,h]*2
		end
	end
	y = Y[:]
	return y
end

function impTTV(G::SparseMatrixCSC{Int64,Int64},H::SparseMatrixCSC{Int64,Int64},x::Vector{Float64},nG::Int,nH::Int)
	X = reshape(x,nG,nH)
	Y = similar(X)
	Y .= 0
	for g = 1:nG
		for h = 1:nH
			g_triangles = triangles(G,g) # no symmetry
			h_triangles = triangles(H,h) # no symmetry
			for g_tri in g_triangles
				j,k = g_tri.w,g_tri.x
				for h_tri in h_triangles
					jp,kp = h_tri.w,h_tri.x
					Y[g,h] += X[j,jp]*X[k,kp]+X[j,kp]*X[k,jp]
				end
			end
			Y[g,h] = Y[g,h]*2
		end
	end
	y = Y[:]
	return y
end

function score_fn(X::Array{Float64,2},A,B)
	ma,mb = edge_list(bipartite_matching(sparse(X)))
	Ap = A[ma,mb]
	Bp = B[mb,mb]
	C = Ap.*Bp
	mytriangles = triangles(A)
	z = collect(mytriangles)
	length(z)
end