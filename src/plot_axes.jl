function PlotAxes.asplotable(x::AuditorySpectrogram,args...;quantize=(100,128),
    kwds...)
  logfreqs = log.(ustrip.(freqs(x)))
  # if it's essentially a range, make it a range
  if std(diff(logfreqs)) < 1e-8
    logfreqs = range(first(logfreqs),last(logfreqs),length=length(logfreqs))
    x = AxisArray(getcontents(x),AxisArrays.axes(x,1),Axis{:logfreq}(logfreqs))
  else
    x = AxisArray(getcontents(x),AxisArrays.axes(x,1),Axis{:logfreq}(logfreqs))
  end

  asplotable(x,args...;quantize=quantize,kwds...)
end

# function PlotAxes.asplotable(x::Cortical,args...;kwds...)
#   logfreqs = Axis{:logfreq}(log.(ustrip.(freqs(x))))
#   x = AxisArray(Array(x),AxisArrays.axes(x)[1:end-1]...,logfreqs)
#   asplotable(x,args...;kwds...)
# end

            
            
