using AxisArrays

export rates, scales, nrates, nscales, default_rates, default_scales,
  cortical, cycoct, co

@dimension Sc "Sc" Scale
@refunit cycoct "cyc/oct" CyclesPerOct Sc false

struct CParams{R,S} <: Params
  aspect::ASParams
  rates::R
  scales::S
  bandonly::Bool

  function CParams(aspect::ASParams,rates::R,scales::S,
                   bandonly::Bool) where {R,S}
    if rates == scales == nothing
      error("You must specify the rates and/or scales.")
    end
    new{R,S}(aspect,
             rates == nothing ? nothing : sort(rates),
             scales == nothing ? nothing :sort(scales),bandonly)
  end
end
const CParamScales{S} = CParams{Void,S}
const CParamRates{R} = CParams{R,Void}
const CParamAll = CParams{R,S} where {R <: AbstractArray,S <: AbstractArray}

struct Cortical{T,N,R,S} <: Result{T,N}
  val::AxisArray{T,N}
  params::CParams{R,S}
end
const CorticalScales{T,N,S} = Cortical{T,N,Void,S}
const CorticalRates{T,N,R} = Cortical{T,N,R,Void}
AxisArrays.AxisArray(x::Cortical) = x.val
Params(x::Cortical) = x.params
similar_helper(::Cortical,data,params) = Cortical(data,params)
resultname(x::Cortical) = "Cortical Rates × Scales"
resultname(x::CorticalRates) = "Cortical Rates"
resultname(x::CorticalScales) = "Cortical Scales"

function showparams(io,params::CParams)
  write(io,"bandonly = $(params.bandonly)\n")
  showparams(io,params.aspect)
end

asrates(x::CParams) = CParams(x.aspect,x.rates,nothing,x.bandonly)
asscales(x::CParams) = CParams(x.aspect,nothing,x.scales,x.bandonly)

cortical_progress(n) = Progress(desc="Cortical Model: ",n)

frame_length(x::Cortical) = frame_length(x.params.aspect)

freqs(x::CParams) = freqs(x.aspect)

rates(x::CParams) = x.rates
rates(x::Result) = rates(AxisArray(x))
rates(x::AxisArray) = axisvalues(axes(x,Axis{:rate}))[1]
nrates(x) = length(rates(x))

scales(x::CParams) = x.scales
scales(x::Result) = scales(AxisArray(x))
scales(x::AxisArray) = axisvalues(axes(x,Axis{:scale}))[1]
nscales(x) = length(scales(x))

Δt(c::CParams) = Δt(c.aspect)
Δf(c::CParams) = Δf(c.aspect)

hastimes(c::Cortical) = HasTimes()

const default_rates = sort([-2.^(1:0.5:5); 2.^(1:0.5:5)]).*Hz
const default_scales = (2.^(-2:0.5:3)).*cycoct

CParams(x::AbstractArray;rates=nothing,scales=nothing,
        bandonly=false,params...) =
  CParams(ASParams(params),rates,scales,bandonly)
CParams(x::AuditorySpectrogram;rates=nothing,scales=nothing,
        bandonly=false) =
  CParams(x.params,rates,scales,bandonly)
function CParams(x::CorticalRates;rates=nothing,scales=nothing,
                 bandonly=false,params...)
  @assert rates == nothing "Already analyzed rates."
  @assert bandonly == x.params.bandonly "`bandonly` value does not match."
  CParams(x.params.aspect,rates,scales,bandonly)
end

function CParams(x::CorticalScales;rates=nothing,scales=nothing,
                 bandonly=false,params...)
  @assert scales == nothing "Already analyzed scales."
  @assert bandonly == x.params.bandonly "`bandonly` value does not match."
  CParams(x.params.aspect,rates,scales,bandonly)
end

const spect_rate = 24
# TODO: implicity convert sound into cortical representation

# cortical responses of rates and scales simultaneously
asHz(x) = x.*Hz
asHz(::Void) = nothing
ascycoct(x) = x.*cycoct
ascycoct(::Void) = nothing
function cortical(y::AbstractArray;progressbar=true,
                  rates_Hz=nothing,rates=asHz(rates_Hz),
                  scales_cycoct=nothing,scales=ascycoct(scales_cycoct),
                  freq_limits_Hz=(),freq_limits=freq_limits_Hz.*Hz,
                  params...)
  if length(freq_limits) > 0
    startHz,stopHz = freq_limits
    y = y[Axis{:freq}(startHz .. stopHz)]
  end

  params = CParams(y;rates=rates,scales=scales,params...)
  cortical(y,params,progressbar)
end

cortical(y::AbstractVector,params::CParams,progressbar=true) =
  cortical(audiospect(y,params.aspect),params,progressbar)
cortical(y::AbstractMatrix,params::CParams,progressbar=true) =
  cortical(audiospect(y,params.aspect),params,progressbar)

####################
# 'identity' functions: converts various arrays that already contain the
# computed cortical representation
function cortical(x::AxisArray{T,4} where T,params::CParamAll,
                  progressbar=true)
  @assert(all(i > 0 for i in indexin(freqs(x),freqs(params))),
          "Frequency channels of array inconsisent with parameters.")
  @assert(all(i > 0 for i in indexin(scales(x),scales(params))),
           "Missing scales in parameters")
  @assert(all(i > 0 for i in indexin(rates(x),rates(params))),
          "Missing rates in parameters")
  Cortical(x,params)
end

function cortical(x::AxisArray{T,3} where T,params::CParamScales,
                  progressbar=true)
  @assert(all(i > 0 for i in indexin(freqs(x),freqs(params))),
          "Frequency channels of array inconsisent with parameters.")
  @assert :rate ∉ axisnames(x) "Unexpectd rate dimension"
  @assert(all(i > 0 for i in indexin(scales(x),scales(params))),
           "Missing scales in parameters")
  Cortical(x,params)
end


function cortical(x::AxisArray{T,3} where T,params::CParamRates,
                  progressbar=true)
  @assert(all(i > 0 for i in indexin(freqs(x),freqs(params))),
          "Frequency channels of array inconsisent with parameters.")
  @assert(all(i > 0 for i in indexin(rates(x),rates(params))),
          "Missing rates in parameters")
  @assert :scale ∉ axisnames(y) "Unexpectd scale dimension"
  Cortical(x,params)
end

function cortical(y::AbstractArray{T,4} where T,params::CParamAll,
                  progressbar=true)
  f = Axis{:freq}(freqs(params.aspect))
  r = Axis{:rate}(params.rates)
  sc = Axis{:scale}(params.scales)
  t = Axis{:time}(times(params.aspect,y))
  Cortical(AxisArray(y,t,r,sc,f),params)
end

function cortical(y::AbstractArray{T,3} where T,params::CParamRates,
                  progressbar=true)
  f = Axis{:freq}(freqs(params.aspect))
  r = Axis{:rate}(params.rates)
  t = Axis{:time}(times(params.aspect,y))
  Cortical(AxisArray(y,t,r,f),params)
end

function cortical(y::AbstractArray{T,3} where T,params::CParamScales,
                  progressbar=true)
  f = Axis{:freq}(freqs(params.aspect))
  sc = Axis{:scale}(params.scales)
  t = Axis{:time}(times(params.aspect,y))
  Cortical(AxisArray(y,t,sc,f),params)
end

####################
# actual cortical computation
function cortical(y::Result,params::CParamAll,progressbar=true)
  progress = progressbar ? cortical_progress(nrates(params)+1) : nothing
  cs = cortical(y,asscales(params),false)
  next!(progress)
  cortical(cs,asrates(params),progressbar,progress)
end

# cortical responses of rates
function cortical(y::Result,params::CParamRates,progressbar=true,
                  progress=progressbar ? cortical_progress(nrates(params)) :
                  nothing)

  if :rate ∈ axisnames(y)
    warning("Rates already analyzed in the input, ",
            "returning this input unmodified.")
  end

  fir = FIRFiltering(y,Axis{:time})

  cr = initcr(y,params)
  for (ri,HR) in enumerate(rate_filters(fir,cr,params))
    cr[Axis{:rate}(ri)] = view(apply(fir,HR),indices(y)...)
    next!(progress)
  end

  Cortical(cr,params)
end
Cortical(cr::AxisArray{T,4} where T,p::CParamRates) =
  Cortical(cr,CParams(p.aspect,p.rates,scales(cr),p.bandonly))

# cortical responses of scales
vecperm(x::AbstractVector,n) = reshape(x,fill(1,n-1)...,:)
function cortical(y::Result,params::CParamScales,progressbar=true,
                  progress=progressbar ? cortical_progress(nscales(params)) :
                  nothing)
  if :scale ∈ axisnames(y)
    warning("Scales already analyzed in the input, returning ",
            "this input unmodified.")
    y
  end
  fir = FIRFiltering(y,Axis{:freq})

  cs = initcr(y,params)
  for (si,HS) in enumerate(scale_filters(fir,cs,params))
    z = apply(fir,conj.(vecperm([HS; zeros(HS)],ndims(y))))
    cs[Axis{:scale}(si)] = view(z,indices(y)...)
    next!(progress)
  end

  Cortical(cs,params)
end
Cortical(cr::AxisArray{T,4} where T,p::CParamScales) =
  Cortical(cr,CParams(p.aspect,rates(cr),p.scales,p.bandonly))

# inverse of cortical rates and scales
function audiospect(cr::Cortical;norm=0.9,progressbar=true)
  @assert(rates(cr) == rates(cr.params),
          "Missing rates, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")
  @assert(scales(cr) == scales(cr.params),
          "Missing scales, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")

  z_cum = FFTCum(cr)

  progress = progressbar ? cortical_progress(nrates(cr)*nscales(cr)) : nothing
  for (ri,HR) in enumerate(rate_filters(z_cum,cr,use_conj=true))
    for (si,HS) in enumerate(scale_filters(z_cum,cr))
      addfft!(z_cum,cr[:,ri,si,:],HR.*[HS; zeros(HS)]')
      next!(progress)
    end
  end

  t = axes(AxisArray(cr),Axis{:time})
  f = axes(AxisArray(cr),Axis{:freq})
  AuditorySpectrogram(AxisArray(normalize!(z_cum,cr,norm),t,f),
                      Params(cr).aspect)
end

# inverse of scales
function audiospect(cr::CorticalScales;norm=0.9,progressbar=true)
  @assert(scales(cr) == scales(cr.params),
          "Missing scales, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")

  z_cum = FFTCum(cr)

  progress = progressbar ? cortical_progress(nscales(cr)) : nothing
  for (si,HS) in enumerate(scale_filters(z_cum,cr))
    addfft!(z_cum,cr[:,si,:],[HS; zeros(HS)]')
    next!(progress)
  end
  t = axes(AxisArray(cr),Axis{:time})
  f = axes(AxisArray(cr),Axis{:freq})

  AuditorySpectrogram(AxisArray(normalize!(z_cum,cr,norm),t,f),
                      cr.params.aspect)
end

# inverse of rates
function audiospect(cr::CorticalRates;norm=0.9,progressbar=true)
  @assert(rates(cr) == rates(cr.params),
          "Missing rates, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")
  z_cum = FFTCum(cr)

  progress = progressbar ? cortical_progress(nrates(cr)) : nothing
  for (ri,HR) in enumerate(rate_filters(z_cum,cr,use_conj=true))
    addfft!(z_cum,cr[:,ri,:],HR)
    next!(progress)
  end
  t = axes(AxisArray(cr),Axis{:time})
  f = axes(AxisArray(cr),Axis{:freq})

  AuditorySpectrogram(AxisArray(normalize!(z_cum,cr,norm),t,f),cr.params.aspect)
end

################################################################################
# private helper functions

function find_fft_dims(y)
  @assert axisdim(y,Axis{:freq}()) == ndims(y)
  @assert axisdim(y,Axis{:time}()) == 1
  find_fft_dims(size(y))
end
find_fft_dims(y::NTuple{N,Int}) where {N} =
  (nextprod([2,3,5],y[1]),y[2:end-1]...,nextprod([2,3,5],y[end]))

struct FIRFiltering{T,N}
  Y::Array{T,N}
  plan
end

function FIRFiltering(y,axis)
  dims = map(axes(y)) do ax
    if axes(y,axis) == ax
      2nextprod([2,3,5],length(ax))
    else
      length(ax)
    end
  end

  along = axisdim(y,axis)
  Y = fft(pad(y,dims),along)
  FIRFiltering(Y,plan_ifft(Y,along))
end
apply(fir::FIRFiltering,H) = fir.plan * (fir.Y .* H)
Base.size(x::FIRFiltering,i...) = size(x.Y,i...)
Base.ndims(x::FIRFiltering) = ndims(x.Y)

function initcr(y,params::CParamRates)
  r = Axis{:rate}(params.rates)
  ax = axes(y)
  newax = ax[1],r,ax[2:end]...

  AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
end

function initcr(y,params::CParamScales)
  s = Axis{:scale}(params.scales)
  ax = axes(y)
  newax = ax[1:end-1]...,s,ax[end]

  AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
end

# TODO: do this for rates as well
reshape_for(v::Array{T,3},cr::AxisArray{T,3}) where T = v
reshape_for(v::Array{T,4},cr::AxisArray{T,4}) where T = v
reshape_for(v::Array{T,3},cr::AxisArray{T,4}) where T =
    reshape(v,ntimes(cr),1,nfreqs(cr))

# keeps track of cumulative sum of FIR filters
# in frequency-space so we can readily normalize the result.
struct FFTCum{T}
  z::Array{Complex{T},2}
  z_cum::Array{Complex{T},2}
  h_cum::Array{T,2}
  plan
end

function FFTCum(cr::Cortical)
  dims = find_fft_dims(size(cr,1,ndims(cr)))
  mult = 1 .+ (Params(cr).rates != nothing,Params(cr).scales != nothing)
  z = zeros(eltype(cr),dims .* mult)

  FFTCum(z,zeros(z),zeros(real(eltype(z)),size(z)...),plan_fft(z))
end

Base.size(x::FFTCum,i...) = size(x.z_cum,i...)
Base.ndims(x::FFTCum) = ndims(x.z_cum)

function addfft!(x::FFTCum,cr,h)
  x.z[1:ntimes(cr),1:nfreqs(cr)] = cr
  Z = x.plan * x.z
  x.h_cum .+= abs2.(h)
  x.z_cum .+= h .* Z
  x
end

function Base.normalize!(x::FFTCum,cr,norm)
  x.h_cum[:,1] .*= 2
  old_sum = sum(x.h_cum[:,nfreqs(cr)])
  x.h_cum .= norm.*x.h_cum + (1 .- norm).*maximum(x.h_cum)
  x.h_cum .*= old_sum ./ sum(x.h_cum[:,nfreqs(cr)])
  x.z_cum ./= x.h_cum

  spectc = (x.plan \ x.z_cum)[1:ntimes(cr),1:nfreqs(cr)]
  max.(real.(2.*spectc),0)
end

pad(x,lens) = pad(x,lens...)
function pad(x,lens::T...) where T <: Number
  @assert all(size(x) .<= lens)
  y = zeros(eltype(x),lens)
  y[indices(x)...] = x
  y
end

# transforms a bandpass frequency response into either a high or low pass
# response (or leaves it untouched)
function askind(H,len,maxi,kind,nonorm)
  if kind == :band
    H
  else
    old_sum = sum(H)
    if kind == :low
      H[1:maxi-1] = 1
    elseif kind == :high
      H[maxi+1:len] = 1
    else
      error("Unexpected filter kind '$kind'.")
    end
    if !nonorm
      H .= H ./ sum(H) .* old_sum
    end

    H
  end
end

function scale_filters(Y,x,params=Params(x))
  N_f = size(Y,ndims(Y)) >> 1
  smin,smax = extrema(scales(x))
  map(scales(x)) do scale
	  scale_filter(ustrip(uconvert(cycoct,scale)), N_f, spect_rate,
                 params.bandonly ? :band :
                 scale == smin ? :low : scale < smax ? :band : :high)
  end
end

# create the frequency-scale filter (filter along spectral axis)
function scale_filter(scale,len,ts,kind)
  # TODO: ./2 should be outside parens, not change for now since this just
  # changes the magnitude of scales and I want to check that all the other
  # changes I've made make sense before fixing this
  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)

  askind(H,len,indmax(H),kind,false)
end

function rate_filters(Y,x,params=Params(x);use_conj=false)
  N_t = size(Y,1) >> 1
  rmin,rmax = extrema(abs.(rates(x)))

  map(rates(x)) do rate
    rate_filter(ustrip(uconvert(Hz,rate)), N_t, Δt(params.aspect),
                params.bandonly ? :band :
                abs(rate) == rmin ? :low :
                abs(rate) < rmax ? :band : :high,use_conj)
  end
end

# create the temporal-rate filter (filter along temporal axis)
function rate_filter(rate,len,Δt,kind,use_conj=false,return_partial=false)
  t = (0:len-1)*ustrip(uconvert(s,Δt))*abs(rate)
  h = @. sin(2π*t) * t^2 * exp(-3.5t)
  h .-= mean(h)

  H0 = view(fft(pad(h,2len)),1:len)
  A = angle.(H0)
  H = abs.(H0)

  maxH,maxi = findmax(H)
  H ./= maxH
  HR = askind(H,len,maxi,kind,true) .* exp.(A*im)

  if use_conj
    HR = conj.(HR)
  end

  if rate >= 0
    HR = pad(HR,2length(HR))
	else
    HR = pad(HR,2length(HR))
		HR[2:end] .= conj.(reverse(HR[2:end]))
		HR[len+1] = abs(HR[len+2])
	end

  if return_partial
    HR,h
  else
    HR
  end
end
