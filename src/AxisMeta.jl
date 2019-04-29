
struct AxisMeta{Axs <: NamedTuple}
    axes::Axs # map from :axis -> metadata
end

# Do I need this function
function validate(meta::AxisMeta,data::MetaArray{<:AxisArray})
 # TODO: make sure the right number of axes are present
end

# TODO: some sort of construction function
# TODO: some sort of accessor function, e.g. by 
# TODO: someway to add a new axis to an existing description    
    