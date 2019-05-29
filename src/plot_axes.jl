function PlotAxes.asplotable(x::AuditorySpectrogram,args...;quantize=(100,128),
    kwds...)
  logfreqs = log.(ustrip.(frequencies(x)))
  # if it's essentially a range, make it a range
  if std(diff(logfreqs)) < 1e-8
    logfreqs = range(first(logfreqs),last(logfreqs),length=length(logfreqs))
    x = AxisArray(getcontents(x),AxisArrays.axes(x,1),Axis{:logfreq}(logfreqs))
  else
    x = AxisArray(getcontents(x),AxisArrays.axes(x,1),Axis{:logfreq}(logfreqs))
  end

  asplotable(x,args...;quantize=quantize,kwds...)
end

default_quantize(x) = (100,)
default_quantize(x,y) = (100,128,)
default_quantize(x,y,args...) = (100,128,fill(10,length(args))...)

function PlotAxes.asplotable(x::MetaAxisArray;kwds...) 
  asplotable(x,axisnames(x)...;kwds...)
end
function PlotAxes.asplotable(x::MetaAxisArray,ax1,axes...;
  quantize=default_quantize(ax1,axes...),kwds...)
  # TODO: handle time or freq not being present

  logfreqs = log.(ustrip.(frequencies(x)))
  otheraxes = AxisArrays.axes(x)[1:end-1]
  others = axisnames(x)[2:end-1]
  if std(diff(logfreqs)) < 1e-8
    logfreqs = range(first(logfreqs),last(logfreqs),length=length(logfreqs))
    x = AxisArray(getcontents(x),otheraxes...,Axis{:logfreq}(logfreqs))
  else
    x = AxisArray(getcontents(x),otheraxes...,Axis{:logfreq}(logfreqs))
  end
  asplotable(x,:time,:logfreq,others...;quantize=quantize,kwds...)
end

            
            
