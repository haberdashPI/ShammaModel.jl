using AxisArrays
using MetaArrays
using DSP
using JLD2
using FileIO
using SampledSignals
using Statistics

import DSP.Filters.freqs

export freqs, times, nfreqs, ntimes, delta_t, delta_f, Δt, Δf, frame_length,
  audiospect, freq_ticks, duration, hastimes, HasTimes, HasNoTimes,
  timedim, describe_axes

########################################
# cochlear filters

(filter::Filter)(x) = DSP.filt(filter.B,filter.A,x)

########################################
# auditory spectrogram data type and its parameters
const fixed_fs = 8000

struct ASParams{T}
  Δt::typeof(1.0s)
  decay_tc::Float64
  nonlinear::Float64
  octave_shift::Float64
  freq_step::Int
  fs::typeof(1.0Hz)
  cochlear::T
end

const WithAxes{Ax} = AxisArray{<:Any,<:Any,<:Any,Ax}
const AuditorySpectrogram = 
  MetaArray{<:WithAxes{<:Tuple{Axis{:time},Axis{:freq}}}, <:ASParams}
const ASParamLike = Union{ASParams,AuditorySpectrogram}

function describe_axes(io::IO,x)
  for ax in AxisArrays.axes(x)
    println(io,string(AxisArrays.axisname(ax))," ",
            string(round(ustrip(ax.val[1]),digits=2)),
            " - ",string(round(ustrip(ax.val[end]),digits=2)),
            string(unit(ax.val[1])))
  end
end

function Base.show(io::IO,::MIME"text/plain",x::AuditorySpectrogram)
  if !get(io, :compact, false)
    println(io,"Auditory Spectrogram")
    describe_axes(io,x)
  else
    println(io,string(duration(x))," AuditorySpectrogram")
  end
end

resultname(x::AuditorySpectrogram) = "Auditory Spectrogram"

num_base_freqs(as::ASParamLike) = length(as.cochlear.filters)-1
nfreqs(as::ASParamLike) = floor(Int,(num_base_freqs(as))/as.freq_step)
nfreqs(x) = length(freqs(x))

function freqs(as::ASParamLike)
  result = (440.0Hz * 2.0.^(((1:num_base_freqs(as)).-31)./24 .+ as.octave_shift))
  result[1:as.freq_step:end]
end
freqs(as::MetaUnion{AxisArray}) = 
  axisvalues(AxisArrays.axes(as,Axis{:freq}))[1]

ntimes(x) = length(times(x))
times(as::MetaUnion{AxisArray}) = axisvalues(AxisArrays.axes(as,Axis{:time}))[1]
times(p::ASParamLike,x::AbstractArray) = (Base.axes(x,1) .- 1) .* Δt(p)

struct HasTimes end
struct HasNoTimes end
hastimes(x::AuditorySpectrogram) = HasTimes()
hastimes(x::MetaUnion{AxisArray}) = :time ∈ axisnames(x) ? HasTimes() : HasNoTimes()
timedim(x::MetaUnion{AxisArray}) = axisdim(x,Axis{:time})

duration(x::AbstractArray) = last(times(x)) - first(times(x))
duration(x::SampleBuf) = nframes(x) / samplerate(x)

delta_t(x) = Δt(x)
delta_f(x) = Δf(x)

frame_length(params::ASParamLike) = floor(Int,Δt(params) * params.fs)
Δt(params::ASParamLike) = params.Δt
Δf(params::ASParamLike) = 2^(1/24)
usamplerate(params::ASParamLike) = uconvert(Hz,params.fs)
SampledSignals.samplerate(params::ASParamLike) = ustrip(uconvert(Hz,params.fs))

Δt(x) = (ts = times(x); ts[2] - ts[1])
usamplerate(x::SampleBuf) = samplerate(x)*Hz

default_sr(x) = 8000.0Hz
default_sr(x::SampleBuf) = usamplerate(x)

function ASParams(x;fs=default_sr(x),delta_t_ms=10,delta_t=delta_t_ms*ms,
                  freq_step=1,Δt=delta_t,decay_tc=8,nonlinear=-2,
                  octave_shift=-1)
  @assert fs == fixed_fs*Hz "The only sample rate supported is $(fixed_fs)Hz"

  p = ASParams(uconvert(s,float(Δt)),Float64(decay_tc),Float64(nonlinear),
               Float64(octave_shift),Int(freq_step),uconvert(Hz,float(fs)),
               cochlear[])

  recommended_length = 2^(4+p.octave_shift)
  if frame_length(p) < recommended_length
    warn("It's recommended that you have a frame length of at least,"*
         " $recommended_length samples, but you have $(frame_length(p)).")
  end
  p
end

########################################
# auditory spectrogram interface
audiospect(x::AbstractArray;progressbar=true,params...) =
  audiospect(x,ASParams(x;params...),progressbar)

####################
# 'identity' conversions (something that's already basically a spectrogram)
function audiospect(x::AbstractArray,params::ASParams,progressbar=true)
  if ndims(x) <= 2 && size(x,2) <= 2
    # the array probably represents a sound
    @warn "Assuming sample rate of input is $(fixed_fs)."
    audiospect(SampleBuf(x,samplerate(params)))
  else
    # the array probably represents a spectrogram
    @assert(ndims(x) == 2,"Input to audiospect must be a sound or a "*
            "previously created spectrogram")
    @assert(size(x,2) == nfreqs(params),
            "Presumed input to be a spectrogram but the "*
            "number of columns do not match the number of frequency channels.")

    f = Axis{:freq}(freqs(params))
    t = Axis{:time}(times(params,x))
    MetaArray(params,AxisArray(x,t,f))
  end
end

function audiospect(x::AxisArray{T,2} where T,params::ASParams,progressbar=true)
  @assert(nfreqs(x) == nfreqs(params),
          "Frequency channels of array and parameters do not match")

  MetaArray(params,x)
end

function audiospect(x::AuditorySpectrogram,params::ASParams,progressbar=true)
  @assert(MetaArrays.getmeta(x) == params,
          "Parameters of spectrogram and input parameters do not match")
  x
end

####################
# the actual computation of a spectrogram
function audiospect(x::SampleBuf,params::ASParams,progressbar=true)
  if usamplerate(x) != fixed_fs*Hz
    error("Unsupported sample rate. Must be $(fixed_fs).")
  end

  frame_len  = frame_length(params)
  N = ceil(Int,length(x) / frame_len) # of frames
  x_ = if length(x) < N*frame_len
    [x; fill(zero(T),N*frame_len - length(x))]
  else
    x
  end
  audiospect_helper(x_,N,params,progressbar)
end

function audiospect_helper(x::AbstractVector{T}, N, params::ASParamLike,
                           progressbar=true, internal_call=false) where {T}
  M = length(params.cochlear.filters)
  Y = fill(zero(T),N, M-1)
  Y_haircell = !internal_call ? nothing : fill(zero(T),length(x),M-1)

  last_haircell = x |>
    params.cochlear.filters[M] |>
    ion_channels(params) |>
    haircell_membrane(params)

  progress = progressbar ? Progress(desc="Auditory Spectrogram: ",M-1) : nothing
  for ch = (M-1):-1:1
    # initial haircell transduction
    y,last_haircell = x |> 
      params.cochlear.filters[ch] |>
      ion_channels(params) |>
      haircell_membrane(params) |>
      lateral_inhibition(last_haircell)

    # recitfication and temporal integration
    Y[:,ch] = y |> rectify |> temporal_integration(params,N)

    # save the intermediate result y if this is an internal call
    if internal_call; Y_haircell[:,ch] = y end
    next!(progress)
  end

  if internal_call
    Y,Y_haircell
  else
    f = Axis{:freq}(freqs(params))
    t = Axis{:time}(times(params,Y))

    # MetaData(params,AxisArray...)
    MetaArray(params,AxisArray(Y[:,1:params.freq_step:end],t,f))
  end
end

########################################
# inverse of auditory spectorgram
function SampledSignals.SampleBuf(y_in::AuditorySpectrogram;
                                  max_iterations=typemax(Int),
                                  target_error=0.05,progressbar=true)
  @assert(max_iterations < typemax(Int) || target_error < Inf,
          "No stopping criterion specified (max_iterations or target_error).")
  M = length(y_in.cochlear.filters)

  # expand y to include all frequencies
  y = zeros(eltype(y_in),size(y_in,1),M-1)
  f_ixs = minimum(freqs(y_in)) .<= freqs(y_in) .<= maximum(freqs(y_in))
  # NOTE: this does not currently support frequency decimation
  # because I don't yet need it
  @assert sum(f_ixs) == size(y_in,2) "Unxpected frequency resolution."
  y[:,f_ixs] = y_in

  # generate initial guess
  x = inv_guess(y_in,y)
  N = ceil(Int,length(x) / frame_length(y_in)) # of frames

  # iteration setup
  min_err = Inf
  min_x = x

  prog = if !progressbar
    nothing
  elseif max_iterations < typemax(Int)
    Progress(max_iterations,"Inverting Spectrogram: ")
  else
    ProgressThresh(target_error,"Inverting Spectrogram: ")
  end

  for iteration in 1:max_iterations
    if min_err < target_error; break end
    y_in.nonlinear == 0 && standardize!(x)

    ŷ,ŷ_haircell = audiospect_helper(x,N,y_in,false,true)
    x = match_x(y_in,x,y,ŷ,ŷ_haircell)

    err = relative_error!(y,ŷ)

    if err < min_err
      min_x = x
      min_err = err
    elseif err-1 > min_err
      randomized_reset!(x)
    end

    x .*= 1.01

    if progressbar
      if max_iterations < typemax(Int)
        ProgressMeter.next!(prog;showvalues =
                            [(:error,string(100round(min_err,digits=4),"%"))])
      else
        ProgressMeter.update!(prog,min_err)
      end
    end
  end
  ProgressMeter.finish!(prog)

  SampleBuf(min_x,samplerate(y_in))
end

########################################
# visualization utility
function freq_ticks(as)
  a = minimum(freqs(as))
  b = maximum(freqs(as))
  step = 0.25

  helper(step) = round.(filter(f -> a <= f*Hz <= b,1000*2.0.^(-3:step:2)),
                        digits=-1)
  fbreaks = helper(step)
  while length(fbreaks) > 7
    fbreaks = helper(step *= 2)
  end

  fs = ustrip.(uconvert.(Hz,freqs(as)))

  findices = mapslices(abs.(fbreaks .- fs'),dims=2) do row
    _, i = findmin(row)
    i
  end

  fbreaks,findices
end


################################################################################
# private helper functions

function ion_channels(params::ASParamLike)
  if params.nonlinear > 0
    x -> 1.0 ./ (1.0 .+ exp.(.-x./params.nonlinear))
  elseif params.nonlinear == 0
    x -> Float64.(x .> 0.0)
  elseif params.nonlinear == -1
    x -> max.(x,0.0)
  elseif params.nonlinear == -2
    identity
    # TODO: implement halfregu
  else
    error("Non linear factor of $fac not supported")
  end
end

function haircell_membrane(params)
  if params.nonlinear  != -2
    β = exp(-1/((1//2)*2^(4+params.octave_shift)))
    y -> filt([1.0],[1.0; -β],y)
  else
    identity
  end
end

rectify(x) = max.(x,0)

function temporal_integration(params::ASParamLike,N)
  frame_len = frame_length(params)
  if !iszero(params.decay_tc)
    α = exp(-1/(params.decay_tc*2^(4+params.octave_shift)))
    function(x)
      y = filt([1.0],[1.0; -α],x)
      @views y[frame_len*(1:N)]
    end
  else    # short-term average
    frame_len == 1 ? identity : y -> mean(reshape(y, frame_len, N),1)'
  end
end

function lateral_inhibition(last_haircell)
  function(y)
    new_y = y .- last_haircell
    new_y,y
  end
end

function inv_guess(params::ASParamLike,y::AbstractMatrix)
  # the initial guess only uses the first 48 channels
  f = ustrip.(uconvert.(Hz,freqs(params)))[1:48]
  steps = 1:frame_length(params)*size(y,1)
  indices = ceil.(Int,steps / frame_length(params))
  t = steps ./ ustrip(uconvert(Hz,params.fs))

  x = dropdims(standardize!(sum(cos.(2π.*f'.*t) .* view(y,indices,1:48),dims=2)),
               dims=2)
end

function match_x(params::ASParamLike,x,y,ŷ,ŷ_haircell)
  M = length(params.cochlear.filters)
  steps = 1:frame_length(params)*size(y,1)
  indices = ceil.(Int,steps / frame_length(params))

  ratios = map(ŷ,y) do ŷ,y
    !iszero(ŷ) ? y ./ ŷ : !iszero(y) ? 2 : 1
  end

  x .= 0
  ch_norm = params.cochlear.norm
  for ch in 1:M-1
    if params.nonlinear == -2
      y1 = ŷ_haircell[:,ch].*view(ratios,indices,ch)
    else
      y1 = ŷ_haircell[:,ch]
      posi = find(y1 .>= 0)
      y1[posi] .*= view(ratios,indices[posi],ch)

      negi = find(y1 .< 0)
      y1[negi] .*= maximum(y1[posi]) / -minimum(y1[negi])
    end

    x .+= reverse(params.cochlear.filters[ch](reverse(y1))) ./ ch_norm
  end

  x
end

function randomized_reset!(x)
  x .= sign.(x) .+ rand(size(x))
  x .-= mean(x)
  x ./= std(x)
  x
end

function relative_error!(y,ŷ)
  y2sum = sum(y.^2)
  ymean = mean(y)

  ŷ .*= ymean/mean(ŷ)
  sum((ŷ .- y).^2) ./ y2sum
end

function standardize!(x)
  x .-= mean(x)
  x ./= std(x)
end
