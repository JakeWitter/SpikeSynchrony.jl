#=
Julia implementation of the van Rossum distance, from
Houghton & Kreuz (2012), [On the efficient calculation of van Rossum distances](https://doi.org/10.3109/0954898X.2012.673048)

Implementation by Geoge Datseris.
=#

export vanRossum, vanRossum_fast

"""
    vanRossum(u, v, τ)
Calculate the van Rossum distance [1] between two spike trains `u, v` with
time parameter `τ` using the exact analytic expression, equation (5)
of [2]. Notice that our normalization is *twice* as much as the original
paper (1.0 instead of 0.5).

`u, v` are expected to be 1D vectors that only list the times of the spikes.

[1] : van Rossum (2001), [A novel spike distance](https://doi.org/10.1162/089976601300014321)

[2] : Houghton & Kreuz (2012), [On the efficient calculation of van Rossum distances](https://doi.org/10.3109/0954898X.2012.673048)
"""
vanRossum(u, v, τ) = sqrt(vR_ds(u, u, τ) + vR_ds(v, v, τ) - 2vR_ds(u, v, τ))
function vR_ds(u, v, τ) # vanRossum double sum
    s = 0.0
    @inbounds @fastmath for i in 1:length(u), j in 1:length(v)
        s += exp(-abs(u[i] - v[j]) / τ)
    end
    return s
end


"""
    vanRossum_fast(u, v, τ) -> d
Calculate the van Rossum distance between spike vectors `u, v`, with
time parameter `τ`. This function uses the fast algorithm proposed by
Houghton and Kreuz [1], (equation (9)), adjusted here by the same 
factor of 2 as `vanRossum`, so its normalization matches that, not
the paper's.

`u, v` are expected to be **sorted** 1D vectors that only list the times of the spikes.
```math
\\begin{aligned}
    [d(\\mathbf{u}, \\mathbf{v} ; \\tau)]^{2} = & n_{1} + n_{2}
    + 2\\sum_{i} \\sum_{j | u_{i} ; u_{j}} e^{-\\left(u_{i}-u_{j}\\right) / \\tau}
    + 2\\sum_{i} \\sum_{j | v_{i}>v_{j}} e^{-\\left(v_{i}-v_{j}\\right) / \\tau} \\\\
    & - 2\\sum_{i} \\sum_{j | u_{i}>v_{j}} e^{-\\left(u_{i}-v_{j}\\right) / \\tau}
    - 2\\sum_{i} \\sum_{j | v_{i}>u_{j}} e^{-\\left(v_{i}-u_{j}\\right) / \\tau}
\\end{aligned}
```

[1] : Houghton & Kreuz (2012), [On the efficient calculation of van Rossum distances](https://doi.org/10.3109/0954898X.2012.673048)
"""
function vanRossum_fast(u::AbstractVector{T1}, v::AbstractVector{T2}, τ) where {T1,T2}
    @assert issorted(u) && issorted(v)

    nu, nv = length(u), length(v)
    mu, mv = markage_vector(u, τ), markage_vector(v, τ)
    # First two factors of Eq. (9):
    Au, Av = sum(mu), sum(mv)
    # Last two factors of Eq. (9):
    Buv, Bvu = cross_term(u, v, mv, τ), cross_term(v, u, mu, τ; strict=false)

    # doubled to match vanRossum's normalization convention
    # max guards against floating point errors
    return sqrt(max(0.0, nu + nv + 2 * Au + 2 * Av - 2 * Buv - 2 * Bvu))
end

"""
    markage_vector(u, τ)
Create the markage vector, Eq. (6) in the Houghton paper.
"""
function markage_vector(u::AbstractVector{T}, τ) where T
    m = zero(u)
    @inbounds @fastmath for i in Iterators.drop(eachindex(u), 1)
        m[i] = (m[i-1] + 1.0) * exp(-(u[i] - u[i-1]) / τ)
    end
    return m
end

"""
    cross_term(u, v, mv, τ)
Calculate Eq. (12) of the Houghton paper.
```math
\\sum_{i} \\sum_{j | u_{i}>v_{j}} e^{-\\left(u_{i}-v_{j}\\right) / \\tau}
=
\\sum_{i} e^{-\\left(u_{i}-v_{J(i)}\\right) / \\tau}\\left[1+m_{J(i)}(\\mathbf{v})\\right]
```
`strict` controls whether an exact coincidence counts in the sum.
`vanRossum_fast` calls `cross_term` once in each direction, once
with `strict=true` and once without, counting coincidences exactly
once.
"""
function cross_term(
    u::AbstractVector{T1}, v::AbstractVector{T2}, mv::AbstractVector{T2}, τ; strict=true
) where {T1,T2}
    s = 0.0
    L = length(v)
    Jprev = 1
    @inbounds @fastmath for i in eachindex(u)
        if strict
            pos = findlast(j -> v[j] < u[i], Jprev:L)
        else
            pos = findlast(j -> v[j] <= u[i], Jprev:L)
        end
        isnothing(pos) && continue # this happens at first and last index
        J = Jprev + pos - 1
        s += exp(-(u[i] - v[J]) / τ) * (1.0 + mv[J])
        Jprev = J
    end
    return s
end
