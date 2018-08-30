using VegaLite
using DataFrames

export vraster

function vraster(z::AbstractMatrix;x=indices(z,1),y=indices(z,2),kwds...)
  Y = ones(x) .* y'
  X = x .* ones(y)'
  df = DataFrame(x = vec(X),y=vec(Y),z = vec(z))
  vraster(df;kwds...)
end

function vraster(df::DataFrame;value=:z,kwds...)
  if any(.!iszero.(imag.(df[value])))
    vraster_plot__complex(df,value=value;kwds...)
  else
    vraster_plot__real(df,value=value;kwds...)
  end
end

function vraster_plot__real(df::DataFrame;value=:z,x=:x,y=:y,maxbins=300,
                            limits=extrema(real.(df[value])),real_suffix=:_real)
  value_r = Symbol(string(value,real_suffix))
  df = copy(df)
  df[value_r] = real.(df[value])

  colormap = if any(df[value_r] .< 0)
    lim = maximum(abs.(limits))
    limits = (-lim,lim)
    "#".*hex.(cmap[:D1])
  else
    "#".*hex.(cmap[:reds])
  end

  nx = length(unique(df[x]))
  ny = length(unique(df[y]))

  # TODO: check on this working when the final, stable VegaLite.jl
  # is out, right now vlfill and vlstroke are undefined, even though
  # they are in the spec.
  df |>
      markrect() |>
    encoding(xquantitative(field=x,vlbin(maxbins=min(nx,maxbins))),
             yquantitative(field=y,vlbin(maxbins=min(ny,maxbins))),
             vlfill(field=value,typ=:quantitative,
                    scale=vlscale(domain=collect(limits),
                                  range=colormap)),
             vlstroke(field=value,typ=:quantitative,
                      scale=vlscale(domain=collect(limits),
                                    range=colormap)))
end

function vraster_plot__complex(df::DataFrame;value=:z,x=:x,y=:y,maxbins=300,
                              phase_suffix=:_phase,abs_suffix=:_abs)
  colormap = "#".*hex.(cmap[:C6])
  df = copy(df)
  z_phase = Symbol(string(value,phase_suffix))
  z_abs = Symbol(string(value,abs_suffix))
  df[z_phase] = angle.(df[value])
  df[z_abs] = abs.(df[value])

  nx = length(unique(df[x]))
  ny = length(unique(df[y]))

  phasescale = vlscale(domain=[-π,π],range=colormap)
  df |>
    markrect() |>
    encoding(xquantitative(field=x,vlbin(maxbins=min(nx,maxbins))),
             yquantitative(field=y,vlbin(maxbins=min(ny,maxbins))),
             opacity=vlopacity(field=z_abs,typ=:quantitative,
                               scale=vlscale(domain=[0,maximum(df[z_abs])])),
             colorquantitative(field=z_abs,scale=phasescale))
end
