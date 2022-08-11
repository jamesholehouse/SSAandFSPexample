# path_to_folder = "/home/jamesholehouse/github/SSAandFSPexample";
# file_path = join([path_to_folder,"/FSP-TM.jl"]);
# include(file_path);
# using .FSPTM
"""
This code coverts an FSP solution to a Julia defined distribution.
Note this could also have been done using a Categorical distribition
with some modifications.
"""


using Distributions, StatsBase

struct FSPd <: DiscreteUnivariateDistribution
    Pn::Vector{Float64}
end

#### Parameters
Distributions.params(d::FSPd) = d.Pn
truncation(d::FSPd) = length(d.Pn)-1

Distributions.pdf(d::FSPd, k::Int) = d.Pn[k+1]
function Distributions.rand(d::FSPd)
    weights = d.Pn
    items = collect(0:1.0:length(d.Pn)-1)
    return sample(items, Weights(weights))
end
Distributions.mean(d::FSPd) = sum([pdf(d,x)*x for x in 0:truncation(d)])
Distributions.var(d::FSPd) = sum([pdf(d,x)*x^2 for x in 0:truncation(d)]) - mean(d)^2
Distributions.logpdf(d::FSPd, k::Int) = log(pdf(d,k))
Distributions.logpdf(d::FSPd, k::Float64) = logpdf(d,floor(Int,k))
Distributions.cdf(d::FSPd, k::Int) = sum([pdf(d,x) for x in 1:k]);
