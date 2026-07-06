using Test, SpikeSynchrony

f = collect(1.0:10.0)

# Distance between identical vectors should be 0, under all q
for q in collect(0:1.5:3.0)
    @test victor_purpura(f, f, q) == 0.0
end

# Distance to empty list should be the length of the non empty list, under all q
for q in collect(0:1.5:3.0)
    @test victor_purpura(f, Float64[], q) == 10
end
@test victor_purpura(Float64[], Float64[], 1.0) == 0.0

# Adding one spike to an identical spike train should have cost 1.0
g = push!(deepcopy(f), 11.0)
@test victor_purpura(f, g, 1.0) == 1.0

s1 = [0.0, 1.0, 2.0]
s2s = Dict([(dt, [0.0, 1.0, 1.0 + dt]) for dt in [0.1, 0.5, 5.0]])

# Distance between trains differing in final spike time reflects the difference
# in that spike time (1+dt vs 2), up to 2.0 for adding and removing a spike
for (dt, s2) in s2s
    @test victor_purpura(s1, s2, 1.0) ≈ min(abs(1.0 - dt), 2.0)
end

# Increasing q increases distance between trains differing in final spike by
# dt, again up to 2.0
δt = 0.5
for q in collect(0.0:1.0:6.0)
    @test victor_purpura(s1, s2s[0.5], q) ≈ min(q * δt, 2.0)
end

sr = sort(rand(3))
# Distance should be symmetric
for q in collect(0:1.5:3)
    @test victor_purpura(s1, sr, q) == victor_purpura(sr, s1, q)
end
