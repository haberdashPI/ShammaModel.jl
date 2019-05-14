using AxisArrays
using FFTW

export rates, scales, nrates, nscales, default_rates, default_scales,
  cortical, cycoct, co

# re-express the spectral and cortical dimensions to have meta data specific to
# an axis and then make it possible to add an axis to an existing array (maybe
# have a flag to allow multiple axes of the same type)

@dimension Sc "Sc" Scale
@refunit cycoct "cyc/oct" CyclesPerOct Sc false

# NOTE: the `low` and `high` fields are used to determine which filters should
# be low- and high-pass, rather than band-pass
struct ScaleAxis
  low::typeof(1cycoct)
  high::typeof(1cycoct)
end

struct RateAxis
  low::typeof(1Hz)
  high::typeof(1Hz)
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

function ratefilt(y::MetaAxisArray, rates; progressbar=true, bandonly=true,
                  progress=progressbar ? cortical_progress(length(rates)) :
                  nothing, rateax=:rate)
  @assert :time in axisnames(y)
  if rateax != :rate && !occursin("rate",string(rateax))
    error("Rate axis name `rateax` must contain the word 'rate'.")
  end
  if rateax in axisnames(y)
    error("Input already has an axis named `$rateax`. Set `rateax` to a",
          " different value to create a second rate axis.")
  end

  fir = FIRFiltering(y,Axis{:time})
  cr = initrates(y,rates,rateax,bandonly)
  for (ri,HR) in enumerate(rate_filters(fir,cr,rateax))
    cr[Axis{rateax}(ri)] = view(apply(fir,HR),Base.axes(y)...)
    next!(progress)
  end

  cr
end

# cortical responses of scales
vecperm(x::AbstractVector,n) = reshape(x,fill(1,n-1)...,:)
function scalefilt(y::MetaAxisArray, scales; progressbar=true, bandonly=true,
                   progress=progressbar ? cortical_progress(length(scales)) :
                   nothing, scaleax=:scale)
  @assert :freq in axisnames(y)
  if scaleax != :scale && !occursin("scale",string(scaleax))
    error("Scale axis name `scaleax` must contain the word 'scale'.")
  end
  if scaleax in axisnames(y)
    error("Input already has an axis named `$scaleax`. Set `scaleax` to a",
          " different value to create a second scale axis.")
  end

  fir = FIRFiltering(y,Axis{:freq})

  cs = initscales(y,scales,scalex,bandonly)
  for (si,HS) in enumerate(scale_filters(fir,cs,scaleax))
    z = apply(fir,conj.(vecperm([HS; zero(HS)],ndims(y))))
    cs[Axis{scaleax}(si)] = view(z,Base.axes(y)...)
    next!(progress)
  end

  cs
end

# TODO: since the normalization step is somewhat inaccurate
# we need some way generic to combine multiple steps

function scaleinv(y::MetaAxisArray;norm=0.9,progressbar=true,scaleax=:scale)
  @assert scaleax in axisnames(y)
 
  z_cum = FFTCum(y)

  progress = progressbar ? cortical_progress(nscales(y)) : nothing
  for (si,HS) in enumerate(scale_filters(z_cum,y,scalex))
    addfft!(z_cum,y[Axis{scaleax}(si)],[HS; zero(HS)]')
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
    reshape(v,ntimes(cr),1,nfreqs(cr))

# keeps track of cumulative sum of FIR filters
# in frequency-space so we can readily normalize the result.
struct FFTCum{T,P,N}
  z::Array{Complex{T},N}
  z_cum::Array{Complex{T},N}
  h_cum::Array{T,N}
  plan::P
end

# TODO: working on generalizing FFTCum to working
# for any number of scale and rate axes
function withoutdim(dims,without)
  dims[1:without-1]...,dims[without+1:end]...
end

function FFTCum(cr::MetaAxisArray,withoutax)
  dims = find_fft_dims((size(cr,1),size(cr,ndims(cr))))
  if occursin("scale",string(withoutax))
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

function scale_filters(Y,x,scaleax)
  N_f = size(Y,ndims(Y)) >> 1
  scale = getproperty(x,scaleax)
  map(scales(x)) do scale
	  scale_filter(ustrip(uconvert(cycoct,scale)), N_f, spect_rate,
                 scale == scale.low ? :low : 
                 scale < scale.high ? :band : :high)
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
  rate = getproperty(x,rateax)

  map(rates(x)) do rate
    rate_filter(ustrip(uconvert(Hz,rate)), N_t, x.time.Δ,
                abs(rate) == rate.low ? :low :
                abs(rate) < rate.high ? :band : :high,use_conj)
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
