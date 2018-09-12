using JLD2
using FileIO
using HDF5

# NOTE: we use jld2 instead of h5 becuase h5 does not easily build in a lot of
# environments (it takes a long time to compile or cannot becuase the
# appropriate compiler isn't available). jld2 builds anywhere that julia does.
cochba = h5open(joinpath("cochba.h5")) do file
  read(file,"/real") + read(file,"/imag")*im
end

jldopen(joinpath("cochba.jld2"),"w") do file
  M = size(cochba,2)
  for ch in 1:M
    p  = floor(Int,real(cochba[1, ch]))
    B  = real(cochba[(0:p).+2, ch])
    A  = imag(cochba[(0:p).+2, ch])
    file["filters/$(string(ch))/A"] = A
    file["filters/$(string(ch))/B"] = B
  end
  file["norm"] = imag(cochba[1, M])
end

