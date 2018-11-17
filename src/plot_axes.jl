function PlotAxes.asplotable(x::AuditorySpectrogram,args...;kwds...)
  logfreqs = Axis{:logfreq}(log.(ustrip.(freqs(x))))
  x = AxisArray(Array(x),AxisArrays.axes(x,1),logfreqs)
  df, axes = asplotable(x,args...;kwds...)
  axes[2].scale = :log

  df, axes
end
            
            
