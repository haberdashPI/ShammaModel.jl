using VegaLite
using DataFrames

function vraster(z::AbstractMatrix;x=Base.axes(z,1),y=Base.axes(z,2),
                 maxbins=500,kwds...)
  Y = ones(x) .* y'
  X = x .* ones(y)'
  df = DataFrame(x = vec(X),y=vec(Y),z = vec(z))
  vraster(df;maxbins=min(maxbins,minimum(size(z))),kwds...)
end

function vraster(df::DataFrame;value=:z,kwds...)
  if any(.!iszero.(imag.(df[value])))
    vraster_plot__complex(df,value=value;kwds...)
  else
    vraster_plot__real(df,value=value;kwds...)
  end
end

function vraster_plot__real(df::DataFrame;value=:z,x=:x,y=:y,maxbins=500,
                            limits=extrema(real.(df[value])),real_suffix=:_real)
  hasnegative = any(x -> real(x) < 0,df[value])
  df |>
    @vlplot(
      :rect,width=300,height=200,
      x={field=x, bin={maxbins=maxbins}, typ="quantitative"},
      y={field=y, bin={maxbins=maxbins}, typ="quantitative"},
      color={field=value, aggregate="mean()", typ="quantitative",
             scale={range=hasnegative ? "diverging" : "heatmap"}},
      config={
        view={stroke="transparent"},
        range={heatmap={scheme="orangered"}}
   })
end

function vraster_plot__complex(df::DataFrame;value=:z,x=:x,y=:y,maxbins=300,
                              phase_suffix=:_phase,abs_suffix=:_abs)
  df = copy(df)
  z_phase = Symbol(string(value,phase_suffix))
  z_abs = Symbol(string(value,abs_suffix))
  df[z_phase] = angle.(df[value])
  df[z_abs] = abs.(df[value])
  max_z_abs = maximum(df[z_abs])

  df |>
    @vlplot(
      :rect,width=300,height=200,
      x={field=x, bin={maxbins=maxbins}, typ="quantitative"},
      y={field=y, bin={maxbins=maxbins}, typ="quantitative"},
      color={field=z_phase, aggregate="mean()", typ="quantitative",
             scale={domain=[-pi,pi]}},
      opacity={aggregate="mean()", field=z_abs, typ="quantitative",
               legend={symbolType="square"},
               scale={range=[0,1]}},
      config={
        view={stroke="transparent"},
        range={heatmap={scheme="sinebow"}}
    })
end

# TODO: spectrogram display not quite working yet
# (handling of log isn't right)

function tovega(as::ShammaModel.AuditorySpectrogram)
  ixs = CartesianIndices(as)
  at(ixs,i) = map(x -> x[i],ixs)

  df = DataFrame(response = vec(as),
                 time = vec(ustrip.(uconvert.(s,times(as)[at(ixs,1)]))),
                 freq = vec(ustrip.(uconvert.(Hz,freqs(as)[at(ixs,2)]))))
  fbreaks,findices = freq_ticks(as)
  p = vraster(df,value=:response,x=:time,y=:freq)
  p.params["encoding"]["x"]["axis"] = Dict("title"=>"Time (s)")
  p.params["encoding"]["y"]["axis"] = Dict("title"=>"Frequency (Hz)")

  p
end
