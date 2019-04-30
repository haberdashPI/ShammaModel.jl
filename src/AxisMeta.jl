
struct AxisMeta{Axs}
    axes::Axs # map from :axis -> metadata
end
Base.getproperty(x::AxisMeta,field::Symbol) = getproperty(x.axes,field)

function validate(meta::AxisMeta,data::AxisArray)
    for field in fieldnames(typeof(Base.getfield(meta,:axes)))
        if field ∉ axisnames(data) 
            error("Extra field `$field` in axis metadata.")
        end
    end

    for axis in axisnames(data)
        if axis ∉ fieldnames(typeof(Base.getfield(meta,:axes)))
            error("Missing field `$axis` in axis metadata.")
        end
    end
end

axismeta(data::AxisArray;keys...) = axismeta(keys.data,data)
axismeta(meta::NamedTuple,data::AxisArray) = axismeta(AxisMeta(meta),data)
function axismeta(meta::AxisMeta,data::AxisArray)
    validate(meta,data)
    MetaArray(meta,data)
end

const MetaAxisArray{A,M} = MetaArray{A,M} where {A<:AxisArray,M<:AxisMeta}
const MetaAxisLike = Union{AxisMeta,MetaAxisArray}
function MetaAxisArray(meta::AxisMeta,data::AxisArray)
    validate(meta,data)
    MetaArray(meta,data)
end

function describe_axes(io::IO,x)
  for ax in AxisArrays.axes(x)
    println(io,string(AxisArrays.axisname(ax))," ",
            string(round(ustrip(ax.val[1]),digits=2)),
            " - ",string(round(ustrip(ax.val[end]),digits=2)),
            string(unit(ax.val[1])))
  end
end

function Base.show(io::IO,::MIME"text/plain",x::MetaAxisArray)
  if !get(io, :compact, false)
    println(io,resultname(x))
  else
    println(io,string(duration(x))," ",resultname(x))
    describe_axes(io,x)
  end
end

# TODO: some sort of accessor function, e.g. by 
# TODO: someway to add a new axis to an existing description    
    