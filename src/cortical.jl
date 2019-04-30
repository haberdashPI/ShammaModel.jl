using AxisArrays
using FFTW

export rates, scales, nrates, nscales, default_rates, default_scales,
  cortical, cycoct, co

# re-express the spectral and cortical dimensions to have meta data specific to
# an axis and then make it possible to add an axis to an existing array (maybe
# have a flag to allow multiple axes of the same type)

@dimension Sc "Sc" Scale
@refunit cycoct "cyc/oct" CyclesPerOct Sc false

struct ScaleAxis
  first::typeof(1cycoct)
  last::typeof(1cycoct)
  bandonly::Bool
end

struct RateAxis
  first::typeof(1Hz)
  last::typeof(1Hz)
  bandonly::Bool
end

cortical_progress(n) = Progress(desc="Cortical Model: ",n)

rates(x::MetaUnion{AxisArray}) =
  axisvalues(AxisArrays.axes(x,Axis{:rate}))[1]
nrates(x) = length(rates(x))

scales(x::MetaUnion{AxisArray}) =
  axisvalues(AxisArrays.axes(x,Axis{:scale}))[1]
nscales(x) = length(scales(x))

const default_rates = sort([-2 .^ (1:0.5:5); 2 .^ (1:0.5:5)]).*Hz
const default_scales = (2 .^ (-2:0.5:3)).*cycoct
const spect_rate = 24

# cortical responses of rates and scales simultaneously
asHz(x) = x*Hz
asHz(x::Quantity) = uconvert(Hz,x)
ascycoct(x) = x*cycoct
ascycoct(x::Quantity) = uconvert(cycoct,x)

function ratec(y::AbstractArray, rates; progressbar=true, bandonly=true,
               progress=progressbar ? cortical_progress(nrates(params)) :
                nothing,
               force=true)
  if any(x->occursin(r"^rate[0-9]+$",x),string.(axisnames(y))
  fir = FIRFiltering(y,Axis{:time})

  cr = initcr(y,params)
  for (ri,HR) in enumerate(rate_filters(fir,cr,params))
    cr[Axis{:rate}(ri)] = view(apply(fir,HR),Base.axes(y)...)
    next!(progress)
  end

  MetaArray(params,cr)
end
MetaArrays.MetaArray(p::CParamRates,cr::AxisArray{T,4} where T) =
  MetaArray(cparams(p.aspect,p.rates,scales(cr),p.bandonly),cr)

# cortical responses of scales
vecperm(x::AbstractVector,n) = reshape(x,fill(1,n-1)...,:)
function cortical(y::Auditory,params::CParamScales,progressbar=true,
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
    z = apply(fir,conj.(vecperm([HS; zero(HS)],ndims(y))))
    cs[Axis{:scale}(si)] = view(z,Base.axes(y)...)
    next!(progress)
  end

  MetaArray(params,cs)
end
MetaArrays.MetaArray(p::CParamScales,cr::AxisArray{T,4} where T) =
  MetaArray(cparams(p.aspect,rates(cr),p.scales,p.bandonly),cr)

# inverse of cortical rates and scales
function audiospect(cr::Cortical;norm=0.9,progressbar=true)
  @assert(rates(cr) == rates(getmeta(cr)),
          "Missing rates, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")
  @assert(scales(cr) == scales(getmeta(cr)),
          "Missing scales, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")

  z_cum = FFTCum(cr)

  progress = progressbar ? cortical_progress(nrates(cr)*nscales(cr)) : nothing
  for (ri,HR) in enumerate(rate_filters(z_cum,cr,use_conj=true))
    for (si,HS) in enumerate(scale_filters(z_cum,cr))
      addfft!(z_cum,cr[:,ri,si,:],HR.*[HS; zero(HS)]')
      next!(progress)
    end
  end

  t = AxisArrays.axes(cr,Axis{:time})
  f = AxisArrays.axes(cr,Axis{:freq})
  audiospect(AxisArray(normalize!(z_cum,cr,norm),t,f), cr.aspect)
end

# inverse of scales
function audiospect(cr::CorticalScales;norm=0.9,progressbar=true)
  @assert(scales(cr) == scales(getmeta(cr)),
          "Missing scales, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")

  z_cum = FFTCum(cr)

  progress = progressbar ? cortical_progress(nscales(cr)) : nothing
  for (si,HS) in enumerate(scale_filters(z_cum,cr))
    addfft!(z_cum,cr[:,si,:],[HS; zero(HS)]')
    next!(progress)
  end
  t = AxisArrays.axes(cr,Axis{:time})
  f = AxisArrays.axes(cr,Axis{:freq})

  MetaArray(cr.aspect,AxisArray(normalize!(z_cum,cr,norm),t,f))
end

# inverse of rates
function audiospect(cr::CorticalRates;norm=0.9,progressbar=true)
  @assert(rates(cr) == rates(getmeta(cr)),
          "Missing rates, this is a slice of the original data."*
          " Slice inversion is currently unsupported.")
  z_cum = FFTCum(cr)

  progress = progressbar ? cortical_progress(nrates(cr)) : nothing
  for (ri,HR) in enumerate(rate_filters(z_cum,cr,use_conj=true))
    addfft!(z_cum,cr[:,ri,:],HR)
    next!(progress)
  end
  t = AxisArrays.axes(cr,Axis{:time})
  f = AxisArrays.axes(cr,Axis{:freq})

  MetaArray(cr.aspect,AxisArray(normalize!(z_cum,cr,norm),t,f))
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
  dims = map(AxisArrays.axes(y)) do ax
    if AxisArrays.axes(y,axis) == ax
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
  ax = AxisArrays.axes(y)
  newax = ax[1],r,ax[2:end]...

  AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
end

function initcr(y,params::CParamScales)
  s = Axis{:scale}(params.scales)
  ax = AxisArrays.axes(y)
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
  dims = find_fft_dims((size(cr,1),size(cr,ndims(cr))))
  mult = 1 .+ (cr.rates != nothing,cr.scales != nothing)
  z = zeros(eltype(cr),dims .* mult)

  FFTCum(z,copy(z),zeros(real(eltype(z)),size(z)...),plan_fft(z))
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

function normalize!(x::FFTCum,cr,norm)
  x.h_cum[:,1] .*= 2
  old_sum = sum(x.h_cum[:,nfreqs(cr)])
  x.h_cum .= norm.*x.h_cum .+ (1 .- norm).*maximum(x.h_cum)
  x.h_cum .*= old_sum ./ sum(view(x.h_cum,:,nfreqs(cr)))
  x.z_cum ./= x.h_cum

  spectc = view((x.plan \ x.z_cum),1:ntimes(cr),1:nfreqs(cr))
  max.(real.(2 .* spectc),0)
end

pad(x,lens) = pad(x,lens...)
function pad(x,lens::T...) where T <: Number
  @assert all(size(x) .<= lens)
  y = zeros(eltype(x),lens)
  y[Base.axes(x)...] = x
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
      H[1:maxi-1] .= 1
    elseif kind == :high
      H[maxi+1:len] .= 1
    else
      error("Unexpected filter kind '$kind'.")
    end
    if !nonorm
      H .= H ./ sum(H) .* old_sum
    end

    H
  end
end

function scale_filters(Y,x,params=x)
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
  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)

  askind(H,len,argmax(H),kind,false)
end

function rate_filters(Y,x,params=x;use_conj=false)
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
