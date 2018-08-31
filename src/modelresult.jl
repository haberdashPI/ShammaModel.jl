using AxisArrays

abstract type Params end
abstract type Result{T,N} <: AbstractArray{T,N} end
# Result children should implement 2 methods
# AxisArray(x) = an AxisArray of the data
# Params(x) = a Params object
# optionally you can implement resultname(x) to make printing prettier

Base.size(x::Result) = size(AxisArray(x))
@inline @Base.propagate_inbounds Base.getindex(x::Result,i...) =
  modelwrap(x,getindex(AxisArray(x),i...))
@inline @Base.propagate_inbounds Base.getindex(x::Result,i::Int...) =
  getindex(AxisArray(x),i...)
@inline @Base.propagate_inbounds Base.setindex!(x::Result,v,i...) =
  setindex!(AxisArray(x),v,i...)
@inline @Base.propagate_inbounds Base.setindex!(x::Result{T},v::T,
                                                i::Int...) where T =
  setindex!(AxisArray(x),v,i...)

modelwrap(x::Result{T},val::T) where T = val
modelwrap(x::A,val::AxisArray{T,N}) where {T,N,A <: Result{T,N}} =
  A(val,Params(x))
modelwrap(x::A,val::AxisArray{T}) where {T,A <: Result{T}} = val
modelwrap(x::Result,val::Result) = val
modelwrap(x::Result,val::T) where T =
  error("unexpected result type $T from getindex on a `Result`.")

Base.IndexStyle(x::Result) = IndexStyle(AxisArray(x))

function Base.similar(x::Result,::Type{S},
                      dims::NTuple{N,Int}) where {S,N}
  if size(x) == dims
    data = similar(AxisArray(x),S)
    similar_helper(x,data,x.params)
  else
    similar(AxisArray(x),S,dims)
  end
end

Base.similar(x::Result,ax1::Axis,axs::Axis...) = similar(x,eltype(x),ax1,axs...)
function Base.similar(x::Result,::Type{S},ax1::Axis,axs::Axis...) where S
  if ndims(x) == 1+length(axs)
    data = similar(AxisArray(x),ax1,axs...)
    similar_helper(x,data,x.params)
  else
    similar(AxisArray(x),S,ax1,axs...)
  end
end

function showparams(io,x)
  for field in fieldnames(typeof(x))
    showfield(io,field,getfield(x,field))
  end
end

function showfield(io,field,value)
  write(io,"$field = ")
  show(io,value)
  write(io,"\n")
end

function showfield(io,field,value::Params)
  showparams(io,value)
end

resultname(x::Result) = string(typeof(x))

function Base.show(io::IO,x::Result{T,N}) where {T,N}
  if get(io, :compact, false)
    write(io,"<")
    write(io,resultname(x))
    write(io,": ")
    write(io,"$(N)d $T data")
    write(io,">")
  else
    write(io,resultname(x))
    write(io,":\n")
    showparams(io,Params(x))
    write(io,"-------\n")
    write(io,"$(N)d $T data with axes: \n")
    for ax in AxisArrays.axes(AxisArray(x))
      write(io,"$(AxisArrays.axisname(ax)): "*
            "$(round(ustrip(ax.val[1]),digits=2))"*
            " - $(round(ustrip(ax.val[end]),digits=2)) $(unit(ax.val[1]))\n")
    end
  end
end
Base.show(io::IO,::MIME"text/plain",x::Result) = show(io,x)

AxisArrays.axisdim(x::Result,ax) = axisdim(AxisArray(x),ax)
AxisArrays.axes(x::Result,i::Int...) = AxisArrays.axes(AxisArray(x),i...)
AxisArrays.axes(x::Result,T::Type{<:Axis}...) = AxisArrays.axes(AxisArray(x),T...)
AxisArrays.axes(x::Result) = AxisArrays.axes(AxisArray(x))
AxisArrays.axisnames(x::Result) = axisnames(AxisArray(x))
AxisArrays.axisvalues(x::Result) = axisvalues(AxisArray(x))
