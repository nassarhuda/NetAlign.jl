function read_edgelist_insmat(filename::String,T::Type)
    f = open(filename)
    header = readline(f)
    header = split(header)
    all_lines = readlines(f)
    lines_split = split.(all_lines)
    index1 = getindex.(lines_split,1)
    ri = parse.(Int,index1).+1
    index3 = getindex.(lines_split,3)
    wi = parse.(T,index3)
    lj = zeros(T,parse(Int,header[1]))
    lj[ri] = wi
    return lj
end

function read_netalign_probem(mainloc::String)
    location = join([mainloc,"_A.smat"])
    isfile(location) ? A = Int.(MatrixNetworks.readSMAT(location)) : A =[]

    location = join([mainloc,"_B.smat"])
    isfile(location) ? B = Int.(MatrixNetworks.readSMAT(location)) : B =[]

    location = join([mainloc,"_L.smat"])
    isfile(location) ? L = MatrixNetworks.readSMAT(location) : L =[]

    location = join([mainloc,"_S.smat"])
    isfile(location) ? S = Int.(MatrixNetworks.readSMAT(location)) : S =[]

    location = join([mainloc,"_li.smat"])
    isfile(location) ? li = read_edgelist_insmat(location,Int) : li =[]

    location = join([mainloc,"_lj.smat"])
    isfile(location) ? lj = read_edgelist_insmat(location,Int) : lj =[]

    location = join([mainloc,"_lw.smat"])
    isfile(location) ? w = read_edgelist_insmat(location,Float64) : w =[]

    return S,w,li,lj,A,B,L
end
"""
`load_netalign_problem`
=============
Read the network alignment problem. For a full list of problems, type `netalign_datasets()`

Input
---
- `probname`: name of the problem from the folder data/

Example
----
`S,w,li,lj,A,B,L = load_netalign_problem("example-overlap")`

- A is the adjacency matrix of the first graph    
- B is the adjacency matrix of the second graph   \n
- L is of size nA x nB, L[i,j] = prior known similarity score between node i in A and node j in B\n
- S is the matrix that encodes the possible squares. For every two possible matches from L, S(i,j) is 1 if that pair contributes a square. Size of (S) = nnz(L) x nnz(L).    
- li,lj,w are triplet format of the matrix L (that encodes prior similarity known)    
"""
function load_netalign_problem(probname)
    mainloc = joinpath(dirname(dirname(@__FILE__)),"data",probname)
    mainlocA = joinpath(dirname(dirname(@__FILE__)),"data",join([probname,"_A.smat"]))
    mainlocS = joinpath(dirname(dirname(@__FILE__)),"data",join([probname,"_S.smat"]))

    if isfile(mainlocA) || isfile(mainlocS)
        S,w,li,lj,A,B,L = read_netalign_probem(mainloc)
    else
        systemerror("Problem name not recognized. Type `netalign_datasets()` for a list of problems",true)
    end
    return S,w,li,lj,A,B,L
end

"""
`netalign_datasets`
=======
Returns the existing datasets that could be used for alignment

Example
------
netalign_datasets()
"""
function netalign_datasets()
    return ["example-2","example-overlap","lcsh2wiki-full","lcsh2wiki-small"]
end
