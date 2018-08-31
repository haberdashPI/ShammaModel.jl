# must be run in Julia 0.6: PerceptualColourMaps is an out of date package with
# some great functionality
import PerceptualColourMaps: cmap
import Colors: RGB
import Colors

open("reds.txt","w") do stream
  for col in hex.(Colors.colormap("Reds"))
    println(stream,col)
  end
end

open("diverging_colors.txt","w") do stream
  for col in hex.(RGB.(cmap("D1")))
    println(stream,col)
  end
end

open("circular_colors.txt","w") do stream
  for col in hex.(RGB.(cmap("C6")))
    println(stream,col)
  end
end

