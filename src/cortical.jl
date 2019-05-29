using AxisArrays
using FFTW

export rates, scales, nrates, nscales, default_rates, default_scales,
  cortical, cycoct, co, scalefilter, ratefilter

# re-express the spectral and cortical dimensions to have meta data specific to
# an axis and then make it possible to add an axis to an existing array (maybe
# have a flag to allow multiple axes of the same type)

@dimension Sc "Sc" Scale
@refunit cycoct "cyc/oct" CyclesPerOct Sc false

# NOTE: the `low` and `high` fields are used to determine which filters should
# be low- and high-pass, rather than band-pass
struct ScaleAxis
  low::typeof(1.0cycoct)
  high::typeof(1.0cycoct)
end

struct RateAxis
  low::typeof(1.0Hz)
  high::typeof(1.0Hz)
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

# cortical responses of rates 
ascycoct(x) = x*cycoct
ascycoct(x::Quantity) = uconvert(cycoct,x)

struct TimeRateFilter
  data::Vector{typeof(1.0s)}
  bandonly::Bool
  axis::Symbol
end

function ratefilter(rates=default_rates;bandonly=true,axis=:rate)
  if axis != :rate && !occursin("rate",string(axis))
    error("Rate axis name `$axis` must contain the word 'rate'.")
  end
  TimeRateFilter(rates,bandonly,axis)
end

function DSP.filt(rates::TimeRateFilter,y::MetaAxisArray; progressbar=true, 
                  progress=progressbar ? 
                    cortical_progress(length(rates.data)) : nothing)
  @assert :time in axisnames(y)
  if rates.axis in axisnames(y)
    error("Input already has an axis named `$(rates.axis)`. If you intended ",
          "to add a second rate dimension, change the `axis` keyword argument ",
          "of `ratefilter` to a different value to create a second rate axis.")
  end

  fir = FIRFiltering(y,Axis{:time})
  cr = initrates(y,rates)
  for (ri,HR) in enumerate(rate_filters(fir,cr,rates.axis))
    cr[Axis{rates.axis}(ri)] = view(apply(fir,HR),Base.axes(y)...)
    next!(progress)
  end

  cr
end

# inverse of rates
struct TimeRateFilterInv
  rates::TimeRateFilter
  norm::Float64
end
Base.inv(rates::TimeRateFilter;norm=0.9) = TimeRateFilter(rates,norm)

function DSP.filt(rateinv::TimeRateFilterInv,cr::MetaAxisArray,progressbar=true)
  @assert rateinv.rates.axis in axisnames(cr)
  z_cum = FFTCum(cr,rateinv.rates.axis)

  progress = progressbar ? cortical_progress(nrates(cr)) : nothing
  for (ri,HR) in enumerate(rate_filters(z_cum,cr,use_conj=true))
    addfft!(z_cum,cr[:,ri,:],HR)
    next!(progress)
  end

  MetaAxisArray(removeaxes(getmeta(cr),rateinv.rates.axis),
    normalize!(z_cum,cr,rateinv.norm))
end

# cortical responses of scales
vecperm(x::AbstractVector,n) = reshape(x,fill(1,n-1)...,:)
struct FreqScaleFilter
  data::Vector{typeof(1.0cycoct)}
  bandonly::Bool
  axis::Symbol
end

function scalefilter(scales=default_scales;bandonly=true,axis=:scale)
  if axis != :scale && !occursin("scale",string(axis))
    error("Scale axis name `$axis` must contain the word 'scale'.")
  end
  FreqScaleFilter(scales,bandonly,axis)
end

function DSP.filt(scales::FreqScaleFilter,y::MetaAxisArray; progressbar=true, 
                   progress=progressbar ? 
                     cortical_progress(length(scales.data)) : nothing)
  @assert :freq in axisnames(y)

  if scales.axis in axisnames(y)
    error("Input already has an axis named `$(scales.axis)`. If you intended ",
          "to add a second rate dimension, change the `axis` keyword argument ",
          "of `scalefilter` to a different value to create a second rate axis.")
  end

  fir = FIRFiltering(y,Axis{:freq})

  cs = initscales(y,scales)
  for (si,HS) in enumerate(scale_filters(fir,cs,scales.axis))
    z = apply(fir,conj.(vecperm([HS; zero(HS)],ndims(y))))
    cs[Axis{scales.axis}(si)] = view(z,Base.axes(y)...)
    next!(progress)
  end

  cs
end

# inverse of scales

struct FreqScaleFilterInv
  scales::FreqScaleFilter
  norm::Float64
end
Base.inv(scales::FreqScaleFilter;norm=0.9) = FreqScaleFilterInv(scales,norm)

function DSP.filt(scaleinv::FreqScaleFilterInv,cr::MetaAxisArray,progressbar=true)
  @assert scaleinv.scales.axis in axisnames(cr)
 
  z_cum = FFTCum(cr,scaleinv.scales.axis)

  progress = progressbar ? cortical_progress(nscales(cr)) : nothing
  for (si,HS) in enumerate(scale_filters(z_cum,cr,scaleinv.scales.axis))
    addfft!(z_cum,cr[Axis{scaleinv.scales.axis}(si)],[HS; zero(HS)]')
    next!(progress)
  end
  MetaAxisArray(removeaxes(getmeta(cr),scaleinv.scales.axis),
    normalize!(z_cum,cr,scaleinv.norm))
end

################################################################################
# private helper functions

function find_fft_dims(y)
  @assert axisdim(y,Axis{:freq}) == ndims(y)
  @assert axisdim(y,Axis{:time}) == 1
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

initrates(y,rates::TimeRateFilter) = 
  initrates(y,rates.data,rates.axis,rates.bandonly)
function initrates(y,rates,rateax=:rate,bandonly=false)
  rates = sort(rates)
  r = Axis{rateax}(rates)
  ax = AxisArrays.axes(y)
  newax = ax[1],r,ax[2:end]...

  axar = AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
  rate_axis = bandonly ? 
     RateAxis(-Inf*Hz,Inf*Hz) : RateAxis(first(rates),last(rates))
  axis_meta = addaxes(getmeta(y);Dict(rateax => rate_axis)...)
  MetaAxisArray(axis_meta,axar)
end

initscales(y,scales::FreqScaleFilter) = 
  initscales(y,scales.data,scales.axis,scales.bandonly)
function initscales(y,scales,scaleax=:scale,bandonly=false)
  scales = sort(scales)
  s = Axis{scaleax}(scales)
  ax = AxisArrays.axes(y)
  newax = ax[1:end-1]...,s,ax[end]

  axar = AxisArray(zeros(complex(eltype(y)),length.(newax)...),newax...)
  scale_axis = bandonly ? 
    ScaleAxis(-Inf*cycoct,Inf*cycoct) :
    ScaleAxis(first(scales),last(scales),bandonly)
  axis_meta = addaxes(getmeta(y);Dict(scaleax => scale_axis)...)
  MetaAxisArray(axis_meta,axar)
end

# TODO: do this for rates as well
reshape_for(v::Array{T,3},cr::AxisArray{T,3}) where T = v
reshape_for(v::Array{T,4},cr::AxisArray{T,4}) where T = v
reshape_for(v::Array{T,3},cr::AxisArray{T,4}) where T =
    reshape(v,ntimes(cr),1,nfrequencies(cr))

# keeps track of cumulative sum of FIR filters
# in frequency-space so we can readily normalize the result.
struct FFTCum{T,P,N,Ax}
  z::Array{Complex{T},N}
  z_cum::Array{Complex{T},N}
  h_cum::Array{T,N}
  plan::P
  axes::Ax
end

# TODO: working on generalizing FFTCum to working
# for any number of scale and rate axes
function withoutdim(dims,without)
  dims[1:without-1]...,dims[without+1:end]...
end

function FFTCum(cr::MetaAxisArray,withoutax)
  withoutd = axisdim(cr,Axis{withoutax})
  dims = find_fft_dims(withoutdim(size(cr),withoutd))
  mult = fill(ndims(cr),1)
  mult[1] += startswith(string(withoutax),"scale")
  mult[end] += startswith(string(withoutax),"rate")
  z = zeros(eltype(cr),(dims .* mult)...)

  newaxes = withoutdim(AxisArrays.axes(cr),withoutd)
  FFTCum(z,copy(z),zeros(real(eltype(z)),size(z)...),plan_fft(z),newaxes)
end

Base.size(x::FFTCum,i...) = size(x.z_cum,i...)
Base.ndims(x::FFTCum) = ndims(x.z_cum)

function addfft!(x::FFTCum,cr,h)
  for I in CartesianIndices(size(x.z)[2:end-1])
    x.z[1:ntimes(cr),I,1:nfrequencies(cr)] = cr[1:ntimes(cr),I,1:nfrequencies(cr)]
    Z = x.plan * x.z[:,I,:]
    x.h_cum[:,I,:] .+= abs2.(h)
    x.z_cum[:,I,:] .+= h .* Z
  end
  x
end

function normalize!(x::FFTCum,cr,norm)
  inner_dims = size(x.z)[2:end-1]
  result = similar(cr,real(eltype(cr)),ntimes(cr),inner_dims...,nfrequencies(cr))
  for I in CartesianIndices(size(x.z)[2:end-1])
    x.h_cum[:,I,1] .*= 2
    old_sum = sum(x.h_cum[:,I,nfrequencies(cr)])
    x.h_cum[:,I,:] .= norm.*x.h_cum[:,I,:] .+ (1 .- norm).*maximum(x.h_cum[:,I,:])
    x.h_cum[:,I,:] .*= old_sum ./ sum(view(x.h_cum[:,I,:],:,nfrequencies(cr)))
    x.z_cum[:,I,:] ./= x.h_cum[:,I,:]

    spectc = view((x.plan \ x.z_cum[:,I,:]),1:ntimes(cr),1:nfrequencies(cr))
    result[:,I,:] .= max.(real.(2 .* spectc),0)
  end

  AxisArray(result,x.axes...)
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

function scale_filters(Y,x,scaleax)
  N_f = size(Y,ndims(Y)) >> 1
  scaleparam = getproperty(x,scaleax)
  map(scales(x)) do scale
	  scale_filter(ustrip(uconvert(cycoct,scale)), N_f, spect_rate,
                 scale == scaleparam.low ? :low : 
                 scale < scaleparam.high ? :band : :high)
  end
end

# create the frequency-scale filter (filter along spectral axis)
function scale_filter(scale,len,ts,kind)
  f2 = ((0:len-1)./len.*ts ./ 2 ./ abs(scale)).^2
  H = f2 .* exp.(1 .- f2)

  askind(H,len,argmax(H),kind,false)
end

function rate_filters(Y,x,rateax;use_conj=false)
  N_t = size(Y,1) >> 1
  rateparam = getproperty(x,rateax)

  map(rates(x)) do rate
    rate_filter(ustrip(uconvert(Hz,rate)), N_t, x.time.Δ,
                abs(rate) == rateparam.low ? :low :
                abs(rate) < rateparam.high ? :band : :high,use_conj)
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
