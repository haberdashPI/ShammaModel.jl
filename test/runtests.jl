using Base.Test
using AuditoryModel
using AxisArrays
using Sounds

x = tone(1kHz,1s,rate=8kHz)
X = audiospect(x)
err = 0.05
x_hat = Sound(X,target_error=err)
as_x_hat = audiospect(x_hat)
@testset "Spectrogram" begin
  @test_throws AssertionError audiospect(tone(1kHz,1s))
  @test mean(X[:,0.9kHz ..1.1kHz]) > mean(X[:,1.9kHz .. 2.1kHz])
  @test sum((as_x_hat .* (mean(X) / mean(as_x_hat)) .- X).^2) ./
    sum(X.^2) <= err

  @test length(freqs(X)) > 0
  @test minimum(freqs(X)) > 0Hz
  @test maximum(freqs(X)) > 3kHz
  @test nfreqs(X) == 128
  @test abs(maximum(times(X)) - duration(x)) < 2Δt(X)
  @test Δt(X) ≈ times(X)[2] - times(X)[1]
  @test Δf(X) ≈ freqs(X)[2]/freqs(X)[1]
  @test 3 < length(freq_ticks(X)[1]) < 7
  @test floor(freq_ticks(X)[1]./10) == freq_ticks(X)[1]./10
end

cr = cortical(X,rates=default_rates,scales=default_scales)
X_hat = audiospect(cr)
@testset "Cortical Model" begin
  @test mean(abs,cr[:,:,:,0.9kHz ..1.1kHz]) >
    mean(abs,cr[:,:,:,1.9kHz .. 2.1kHz])
  @test mean((X_hat .- X).^2 ./ sum(X.^2)) < 1e-4
  @test rates(cr) == default_rates
  @test scales(cr) == default_scales
  @test freq_ticks(cr) == freq_ticks(X)
end
