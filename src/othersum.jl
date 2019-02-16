"""
v = othersum(si,s,m) returns an m-by-1 vector where
v[si[i]] = sum(s[1:i])-s[i]
"""  
function othersum(si::Vector{Int64},s::Vector{Float64},m::Int64)
  rowsum = zeros(Float64,m)
  for i = 1:length(si)
    rowsum[si[i]] += s[i]
  end
  os = rowsum[si] .- s
  return os
end