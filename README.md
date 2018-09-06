# ShammaModel

A somewhat minimal implementation of the auditory spectrogram (`audiospect`)
and cortical model (`cortical`) as described in the following paper:

[Chi, T., Ru, P., & Shamma, S. A. (2005). Multiresolution spectrotemporal
analysis of complex sounds. The Journal of the Acoustical Society of America,
118(2), 887–906.](http://doi.org/10.1121/1.1945807)

You can find the original MATLAB implementation of these models
[here](https://isr.umd.edu/Labs/NSL/Software.htm).

# TODO

- Document function interfaces
- Improve plotting (don't require R)
- And a few example uses
- Improve and generify the `ShammaModel.Result` type and make it its own
  package for handling arrays with metadata.

# Status

These functions are quite stable; I'm using them extensively in other projects.
They have not been cleaned up to make them accessible in any way yet.
