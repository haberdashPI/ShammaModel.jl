module ShammaModel
using DataFrames
using Unitful
using Unitful: ms, s, Hz, kHz
using ProgressMeter
using Requires

export ms, s, Hz, kHz

next!(x::Progress) = ProgressMeter.next!(x)
next!(x::Nothing) = nothing

include("modelresult.jl")
include("audiospect.jl")
include("cortical.jl")

const Filter = NamedTuple{(:A, :B),Tuple{Array{Float64,1},Array{Float64,1}}}
const Filters = NamedTuple{(:norm, :filters),Tuple{Float64,Array{Filter}}}
const cochlear = Ref{Filters}()

# include("rplots.jl")
# include("vplots.jl")
const localunits = Unitful.basefactors
const localpromotion = Unitful.promotion
function __init__()
  cochlear[] = load(joinpath(@__DIR__,"..","data","cochba.jld2"),"cochba")
  @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" include("rplots.jl")
  # @require VegaLite include("rplots.jl")

  merge!(Unitful.basefactors, localunits)
  merge!(Unitful.promotion, localpromotion)
  Unitful.register(ShammaModel)
end

end
