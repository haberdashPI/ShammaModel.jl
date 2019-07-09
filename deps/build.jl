using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"
  Pkg.Registry.add("https://github.com/haberdashPI/LabRegistry.jl")
end
