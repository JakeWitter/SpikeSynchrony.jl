#=
Julia implementation of the Victor-Purpura metric from
Victor & Purpura (1997), [Metric-space analysis of spike trains: theory, algorithms and application](https://arxiv.org/abs/q-bio/0309031)
implementation by Jake Witter
=#
export victor_purpura

"""
    victor_purpura(u, v, q) -> d

Calculate the Victor-Purpura distance [1] between spike trains `u, v`, with time
scale parameter `q`. This is computed via the dynamic programming algorithm of
equation (4), applied to spike times rather than intervals, as in the
"Distances based on spike times" section. It uses insertion and deletion costs
of 1.

`u, v` are expected to be 1D sorted vectors listing the times of spikes.

[1] : Victor & Purpura (1997), Metric-space analysis of spike trains: theory,
algorithms, and application
"""
function victor_purpura(u::AbstractVector{T1}, v::AbstractVector{T2}, q::Real) where {T1,T2}
    @assert issorted(u) && issorted(v)

    m = length(u)
    n = length(v)
    G = zeros(m + 1, n + 1)
    G[:, 1] = 0:m
    G[1, :] = 0:n

    for i = 2:m+1
        for j = 2:n+1
            G[i, j] = min(G[i-1, j] + 1,
                G[i, j-1] + 1,
                G[i-1, j-1] + q * abs(u[i-1] - v[j-1]))
        end
    end
    return G[m+1, n+1]
end
