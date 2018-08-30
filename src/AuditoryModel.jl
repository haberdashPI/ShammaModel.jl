module AuditoryModel
using DataFrames
using Unitful
using Unitful: ms, s, Hz, kHz
using ProgressMeter
using Requires

export ms, s, Hz, kHz

next!(x::Progress) = ProgressMeter.next!(x)
next!(x::Void) = nothing

include("modelresult.jl")
include("audiospect.jl")
include("cortical.jl")

@require RCall include(joinpath(@__DIR__,"rplots.jl"))
 # include(joinpath(@__DIR__,"rplots.jl"))
# @require VegaLite include(joinpath(@__DIR__,"rplots.jl"))
# include(joinpath(@__DIR__,"vplots.jl"))

end
