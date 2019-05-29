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
x_hat = filt(AudiospectInv(target_error=err,max_iterations=1000),X)
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

@testset "Data is plottable" begin
  @test size(PlotAxes.asplotable(X)[1],1) == length(X)
end

cr = cortical(X,rates=default_rates,scales=default_scales)
X_hat = audiospect(cr)
@testset "Cortical Model" begin
  @test mean(abs,cr[:,:,:,0.9kHz ..1.1kHz]) >
    mean(abs,cr[:,:,:,1.9kHz .. 2.1kHz])
  @test quantile(vec((X_hat .- X).^2 ./ mean(abs2,X)),0.75) < 5e-3
  @test rates(cr) == default_rates
  @test scales(cr) == default_scales
  @test freq_ticks(cr) == freq_ticks(X)
  @test Array(cortical(cortical(X,scales=default_scales),rates=default_rates)) ≈
    Array(cr)
end

