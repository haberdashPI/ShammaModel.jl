function PlotAxes.asplotable(x::AuditorySpectrogram,args...;kwds...)
  logfreqs = Axis{:logfreq}(log.(ustrip.(freqs(x))))
  x = AxisArray(Array(x),AxisArrays.axes(x,1),logfreqs)
  asplotable(x,args...;kwds...)
end
            
            
