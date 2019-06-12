using Test
using AxisArrays
using SampledSignals
using Statistics
using Unitful
using DSP
using PlotAxes
using ShammaModel

x = SampleBuf(sin.(2π .* 1000 .* range(0,stop=1,length=8000)),8000)
X = filt(audiospect,x)

err = 0.05
x_hat = filt(inv(audiospect,target_error=err,max_iterations=1000),X)
as_x_hat = filt(audiospect,x_hat)

@testset "Spectrogram" begin
  @test_throws ErrorException filt(audiospect,SampleBuf(collect(1:10),4000))
  @test mean(X[:,0.9kHz ..1.1kHz]) > mean(X[:,1.9kHz .. 2.1kHz])
  @test sum((as_x_hat .* (mean(X) / mean(as_x_hat)) .- X).^2) ./
    sum(X.^2) <= err

  @test eltype(filt(audiospect,collect(1:10))) == float(Int)

  @test length(frequencies(X)) > 0
  @test minimum(frequencies(X)) > 0Hz
  @test maximum(frequencies(X)) > 3kHz
  @test nfrequencies(X) == 128
  @test nfrequencies(X[:,1:20]) == 20
  @test abs(maximum(times(X)) - duration(x)*s) < 2Δt(X)
  @test Δt(X) ≈ times(X)[2] - times(X)[1]
  @test Δf(X) ≈ frequencies(X)[2]/frequencies(X)[1]
  @test 3 < length(freq_ticks(X)[1]) < 7
  @test floor.(freq_ticks(X)[1]./10) == freq_ticks(X)[1]./10
end

# TODO: add tests for plotting cortical data as well
@testset "Data is plottable" begin
  @test size(PlotAxes.asplotable(X)[1],1) == length(X)
end

@testset "Single dimension cortical model" begin
  scalef = scalefilter()
  S_cr = filt(scalef,X)
  S_X̂ = filt(inv(scalef),S_cr)

  ratef = ratefilter()
  R_cr = filt(ratef,X)
  R_X̂ = filt(inv(ratef),R_cr)

  @test mean(abs,S_cr[:,:,0.9kHz ..1.1kHz]) >
    mean(abs,S_cr[:,:,1.9kHz .. 2.1kHz])
  @test mean(abs,R_cr[:,:,0.9kHz ..1.1kHz]) >
    mean(abs,R_cr[:,:,1.9kHz .. 2.1kHz])
  @test quantile(vec((S_X̂ .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test quantile(vec((R_X̂ .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test isempty(setdiff(scales(S_cr),default_scales))
  @test isempty(setdiff(rates(R_cr),default_rates))
  @test freq_ticks(S_cr) == freq_ticks(X)
  @test freq_ticks(R_cr) == freq_ticks(X)
end

@testset "Multi dimensional cortical model" begin
  cort = cortical()
  cr = filt(cort,X)
  X̂ = filt(inv(cort),cr)

  # TODO: fix tests
  @test mean(abs,cr[:,:,:,0.9kHz ..1.1kHz]) >
    mean(abs,cr[:,:,:,1.9kHz .. 2.1kHz])
  @test quantile(vec((X̂ .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test isempty(setdiff(scales(cr),default_scales))
  @test isempty(setdiff(rates(cr),default_rates))
  @test freq_ticks(cr) == freq_ticks(X)
end

@testset "Repeated axis cortical model" begin
  cort = cortical()
  cr = filt(cort,X)

  ratef = ratefilter(axis=:rateslow)
  crR = filt(ratef,cr)
  @test ndims(crR) == 5

  ratef = ratefilter()
  @test_throws(ArgumentError("axis name :rate is used more than once"),
    filt(ratef,cr))
end
