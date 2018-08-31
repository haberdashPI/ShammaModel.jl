using JLD2
using FileIO
using HDF5

# NOTE: we use jld2 instead of h5 becuase h5 does not easily build in a lot of
# environments (it takes a long time to compile or cannot becuase the
# appropriate compiler isn't available). jld2 builds anywhere that julia does.
cochba = h5open(joinpath("data","cochba.h5")) do file
  read(file,"/real") + read(file,"/imag")*im
end

M = size(cochba,2)
filters = map(1:M) do ch
  p  = floor(Int,real(cochba[1, ch]))
  B  = real(cochba[(0:p)+2, ch])
  A  = imag(cochba[(0:p)+2, ch])

  (A=A,B=B)
end

ch_norm  = imag(cochba[1, M])
data = (norm=ch_norm,filters=filters)
save(joinpath("data","cochba.jld2"),"cochba",data)

