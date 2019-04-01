# ShammaModel

A somewhat minimal implementation of the auditory spectrogram (`audiospect`)
and cortical model (`cortical`) as described in the following paper:

[Chi, T., Ru, P., & Shamma, S. A. (2005). Multiresolution spectrotemporal
analysis of complex sounds. The Journal of the Acoustical Society of America,
118(2), 887–906.](http://doi.org/10.1121/1.1945807)

You can find the original MATLAB implementation of these models
[here](https://isr.umd.edu/Labs/NSL/Software.htm).

## Status

These functions are quite stable; I'm using them extensively in another project.
I have not put much work into documenting the interface yet, but it is relatively
straighforward to use. Take a look at `test/runtests.jl` for some examples.

## TODO

- Document function interfaces
- And a few example uses, and usage with `PlotAxes` for visualization.
