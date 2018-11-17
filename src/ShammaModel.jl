module ShammaModel
using DataFrames
using Unitful
using Unitful: ms, s, Hz, kHz
using ProgressMeter
using Requires

export ms, s, Hz, kHz

next!(x::Progress) = ProgressMeter.next!(x)
next!(x::Nothing) = nothing

struct Filter
  A::Vector{Float64}
  B::Vector{Float64}
end
struct CochFilters
  filters::Vector{Filter}
  norm::Float64
end
const cochlear = Ref{CochFilters}()

include("audiospect.jl")
include("cortical.jl")

# include("rplots.jl")
# include("vplots.jl")
const localunits = Unitful.basefactors
const localpromotion = Unitful.promotion
function __init__()
  cochlear[] = jldopen(joinpath(@__DIR__,"..","data","cochba.jld2"),"r") do file
    filters = map(keys(file["filters"])) do ch
      A = file["filters/$ch/A"]
      B = file["filters/$ch/B"]
      Filter(A,B)
    end
    CochFilters(filters,file["norm"])
  end

  @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" include("rplots.jl")
  @require PlotAxes="8b6f5f00-d239-11e8-3a24-33314b00f6b0" begin
    using .PlotAxes
    include("plot_axes.jl")
  end
  # @require VegaLite include("rplots.jl")

  merge!(Unitful.basefactors, localunits)
  merge!(Unitful.promotion, localpromotion)
  Unitful.register(ShammaModel)
end

end
