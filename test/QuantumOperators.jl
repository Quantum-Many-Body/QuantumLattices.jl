using LaTeXStrings: latexstring
using LinearAlgebra: dot
using Printf: @sprintf
using QuantumLattices: ZeroAtLeast, ⊗, add!, div!, expand, mul!, id, rank, sub!, update!, value
using QuantumLattices.QuantumOperators
using QuantumLattices.Toolkit: Float, contentnames, isparameterbound, parameternames, parametertype, subscript, superscript

import QuantumLattices: permute
import QuantumLattices.QuantumOperators: script

struct AID{O<:Real, S<:Real} <: OperatorIndex
    orbital::O
    nambu::S
end
@inline Base.adjoint(id::AID) = AID(id.orbital, 3-id.nambu)
@inline script(id::AID, ::Val{:orbital}; kwargs...) = id.orbital
@inline script(id::AID, ::Val{:nambu}; kwargs...) = id.nambu==2 ? "\\dagger" : ""
latexformat(AID, LaTeX{(:nambu,), (:orbital,)}('c'))
function permute(u₁::AID, u₂::AID)
    @assert u₁ ≠ u₂ "permute error: permuted operator units should not be equal to each other."
    if (u₁.nambu == 3-u₂.nambu) && (u₁.orbital == u₂.orbital)
        if u₁.nambu == 2
            return (Operator(1), Operator(1, u₂, u₁))
        else
            return (Operator(-1), Operator(1, u₂, u₁))
        end
    else
        return (Operator(1, u₂, u₁),)
    end
end

@testset "scalartype" begin
    @test scalartype(1) == scalartype(Int) == Int
    @test scalartype(Vector{Int}) == scalartype(Matrix{Int}) == scalartype(Vector{Matrix{Int}}) == Int
    @test scalartype(Tuple{Int}) == scalartype(Tuple{Tuple{Int, Int}}) == Int
end

@testset "ID" begin
    @test promote_type(ZeroAtLeast{AID{Int, Int}, 1}, ZeroAtLeast{AID{Int, Int}}) == ZeroAtLeast{AID{Int, Int}}
    @test promote_type(ZeroAtLeast{AID{Int, Int}}, ZeroAtLeast{AID{Int, Int}, 1}) == ZeroAtLeast{AID{Int, Int}}
    @test promote_type(ZeroAtLeast{AID{Int, Int}}, ZeroAtLeast{AID{Int, Float}}) == ZeroAtLeast{AID{Int, <:Real}}
    @test promote_type(Tuple{}, ZeroAtLeast{AID{Int, Int}}) == ZeroAtLeast{AID{Int, Int}}
    @test promote_type(Tuple{}, ZeroAtLeast{AID{Int, Int}, 2}) == ZeroAtLeast{AID{Int, Int}}

    sid = AID(1, 1)
    @test string(sid) == "AID(1, 1)"
    @test !iszero(sid)
    @test sid'==AID(1, 2)

    cid = (AID(2, 1), AID(Inf, 1))
    @test cid == (AID(2, 1) ⊗ AID(Inf, 1))
    @test ⊗(sid, cid) == ⊗(AID(1, 1), AID(2, 1), AID(Inf, 1))
    @test ⊗(cid, sid) == ⊗(AID(2, 1), AID(Inf, 1), AID(1, 1))
    @test ⊗(cid, cid) == ⊗(AID(2, 1), AID(Inf, 1), AID(2, 1), AID(Inf, 1))
    @test cid == ZeroAtLeast(AID, (2, Inf), (1, 1))
    @test cid|>propertynames == (:orbitals, :nambus)
    @test cid.orbitals == (2, Inf)
    @test cid.nambus == (1, 1)
    @test cid|>eltype == AID{<:Real, Int}
    @test cid|>rank == cid|>typeof|>rank == 2
    @test cid'==(AID(Inf, 2), AID(2, 2))
    @test ishermitian(cid)==false
    @test ishermitian(cid ⊗ cid')
end

@testset "QuantumOperator" begin
    @test contentnames(OperatorProd) == (:value, :id)
    @test parameternames(OperatorProd) == (:value, :id)
    @test isparameterbound(OperatorProd, :value, Any) == false

    @test valtype(OperatorProd) == parametertype(OperatorProd, :value) == parametertype(OperatorProd, 1) == Any
    @test valtype(OperatorProd{Int}) == parametertype(OperatorProd{Int}, :value) == parametertype(OperatorProd{Int}, 1) == Int
    @test idtype(OperatorProd{<:Number}) == parametertype(OperatorProd{<:Number}, :id) == parametertype(OperatorProd{<:Number}, 2) == Tuple
    @test idtype(OperatorProd{<:Number, ZeroAtLeast{AID}}) == parametertype(OperatorProd{<:Number, ZeroAtLeast{AID}}, :id) == parametertype(OperatorProd{<:Number, ZeroAtLeast{AID}}, 2) == ZeroAtLeast{AID}

    @test promote_type(Operator{Int}, Operator) == Operator
    @test promote_type(Operator, Operator{Int}) == Operator
    @test promote_type(Operator{Int}, Operator{Float}) == Operator{Float}
    @test promote_type(Operator{Int, ZeroAtLeast{AID{Int, Int}, 2}}, Operator{Float, ZeroAtLeast{AID{Int, Int}, 2}}) == Operator{Float, ZeroAtLeast{AID{Int, Int}, 2}}
    @test promote_type(Operator{Int, ZeroAtLeast{AID}}, Operator{Float, ZeroAtLeast{AID}}) == Operator{Float, <:ZeroAtLeast{AID}}
    @test promote_type(Operator{Int}, Float) == Operator{Float}

    opt = Operator(2.0, AID(1, 1))
    @test deepcopy(opt)==expand(opt) && isequal(deepcopy(opt), expand(opt))
    @test eltype(opt) == eltype(typeof(opt)) == AID{Int, Int}
    @test collect(opt) == [AID(1, 1)]
    @test opt|>rank == opt|>typeof|>rank == 1
    @test opt|>valtype == opt|>typeof|>valtype == parametertype(opt|>typeof, :value) == Float
    @test opt|>idtype == opt|>typeof|>idtype == parametertype(opt|>typeof, :id) == ZeroAtLeast{AID{Int, Int}, 1}
    @test opt|>scalartype == opt|>typeof|>scalartype == Float
    @test opt|>isequivalenttoscalar == opt|>typeof|>isequivalenttoscalar == false
    @test value(opt) == 2.0
    @test id(opt) == ⊗(AID(1, 1))
    @test replace(opt, 3) == Operator(3, AID(1, 1))
    @test !iszero(opt) && iszero(replace(opt, 0))
    @test isapprox(opt, replace(opt, opt.value+10^-6); atol=10^-5)
    @test length(opt) == 1 && firstindex(opt) == 1 && lastindex(opt) ==1
    @test opt[1]==opt[begin]==opt[end]==AID(1, 1) && opt[1:1]==Operator(1.0, AID(1, 1))
    @test split(opt) == (2.0, AID(1, 1))
    @test opt|>typeof|>one == opt|>one == Operator(1.0)
    @test convert(Operator{Float, Tuple{AID{Int, Int}}}, Operator(2, AID(1, 1))) == Operator(2.0, AID(1, 1))
    @test Operator(2.0) == Operator(2.0)
    @test string(opt) == "Operator(2.0, AID(1, 1))"
    @test opt' == Operator(2.0, AID(1, 2))
    @test ishermitian(opt) == false
    @test ishermitian(Operator(2.0, AID(1, 1), AID(1, 2)))
    @test sequence(opt, Dict(AID(1, 1)=>1, AID(1, 2)=>2)) == (1,)
    @test convert(typeof(opt), AID(1, 1)) == Operator(1.0, AID(1, 1))
    @test convert(Operator{Float, <:ZeroAtLeast{AID}}, AID(1, 1)) == Operator(1.0, AID(1, 1))
    @test convert(Operator{Float, Tuple{}}, 2) == Operator(2.0)
    @test convert(Operator{Float, <:ZeroAtLeast{AID}}, 2) == Operator(2.0)

    opt₁ = Operator(1.0, AID(1, 1))
    opt₂ = Operator(2.0, AID(1, 2))
    opts = Operators(opt₁, opt₂)
    @test opts == Operators{eltype(opts)}(opt₁, opt₂) == OperatorSum(opt₁, opt₂) == OperatorSum((opt₁, opt₂)) == OperatorSum{eltype(opts)}(opt₁, opt₂) == OperatorSum{eltype(opts)}((opt₁, opt₂)) == expand(opts)
    @test eltype(opts) == eltype(typeof(opts)) == Operator{Float, ZeroAtLeast{AID{Int, Int}, 1}}
    @test scalartype(opts) == scalartype(typeof(opts)) == Float
    @test isequivalenttoscalar(opts) == isequivalenttoscalar(typeof(opts)) == false
    @test update!(opts) == opts
    @test collect(opts) == collect(values(opts.contents))
    @test length(opts) == 2
    @test summary(opts) == "Operators"
    @test string(opts) == @sprintf "Operators with 2 Operator\n  %s\n" join(opts, "\n  ")
    @test haskey(opts, id(opt₁)) && haskey(opts, id(opt₂)) && !haskey(opts, ⊗(AID(3, 1)))
    @test opts[1]==opts[begin]==opt₁ && opts[2]==opts[end]==opt₂
    @test opts[1:2] == opts[:] == opts
    @test empty(opts) == empty!(deepcopy(opts)) == zero(opts)
    @test !iszero(opts) && iszero(zero(opts))
    @test !invoke(iszero, Tuple{OperatorSet}, opts) && invoke(iszero, Tuple{OperatorSet}, zero(opts))
    optp₁ = promote_type(typeof(opts), OperatorSum{Operator{Complex{Int}, NTuple{2, AID{Float, Float}}}, NTuple{2, AID{Float, Float}}})
    optp₂ = OperatorSum{Operator{Complex{Float}, <:Tuple{AID, Vararg{AID{Float, Float}}}}, Tuple{AID, Vararg{AID{Float, Float}}}}
    @test optp₁ == optp₂
    @test isapprox(opts, deepcopy(opts)) && isapprox(opts, opts+10^-6; atol=10^-5) && isapprox(opts+10^-6, opts; atol=10^-5)
    @test !isapprox(opts, opts+10^-6) && !isapprox(opts+10^-6, opts)
    @test opts' == Operators(opt₁', opt₂')
    @test !ishermitian(opts) && ishermitian(Operators(opt₁, opt₂, opt₁', opt₂'))
    @test zero(opts) == zero(typeof(opts)) == OperatorSum{eltype(opts)}()
    @test add!(deepcopy(opts))==opts && add!(zero(opts), opts)==opts
    @test sub!(deepcopy(opts))==opts && sub!(deepcopy(opts), opts)==zero(opts)
    @test mul!(deepcopy(opts), 2.0) == opts*2
    @test div!(deepcopy(opts), 2.0) == opts/2

    @test operatortype(AID(1, 1)) == operatortype(AID{Int, Int}) == Operator{Int, Tuple{AID{Int, Int}}}
    @test operatortype(AID) == Operator{Int, <:Tuple{AID}}
    @test operatortype(opt) == operatortype(typeof(opt)) == Operator{Float, Tuple{AID{Int, Int}}}
    @test operatortype(opts) == operatortype(typeof(opts)) == eltype(opts)

    @test zero(AID(1, 1)) == zero(Operator(1, AID(1, 1))) == zero(Operator(1, AID(1, 1))+Operator(1, AID(2, 1)))
    @test conj(AID(1, 1)) == AID(1, 1)
    @test conj(Operator(2im, AID(2, 1), AID(1, 1))) == Operator(-2im, AID(2, 1), AID(1, 1))
    @test conj(Operator(2im, AID(2, 1), AID(1, 1))+Operator(2, AID(1, 1), AID(1, 1))) == Operator(-2im, AID(2, 1), AID(1, 1))+Operator(2, AID(1, 1), AID(1, 1))
    @test dot(Operator(2im, AID(1, 1)), Operator(2im, AID(1, 1), AID(2, 1))) == Operator(4, AID(1, 1), AID(1, 1), AID(2, 1))
    @test dot(Operator(2im, AID(2, 1), AID(1, 1)), 2) == Operator(-4im, AID(2, 1), AID(1, 1))
    @test dot(2im, Operator(2im, AID(2, 1), AID(1, 1))) == Operator(4, AID(2, 1), AID(1, 1))

    @test +opt == opt
    @test -opt == Operator(-2.0, AID(1, 1))
    @test opt*2 == 2*opt == Operator(4.0, AID(1, 1))
    @test opt/2 == Operator(1.0, AID(1, 1))
    @test opt^2 == opt*opt == Operator(4.0, AID(1, 1), AID(1, 1))
    @test opt+1 == 1+opt == Operators(Operator(1), opt)
    @test 1-opt == Operators(Operator(1), -opt)
    @test opt-1 == Operators(opt, Operator(-1))
    @test 2*AID(1, 1) == AID(1, 1)*2 == Operator(2, AID(1, 1))
    @test AID(1, 1)*AID(1, 2) == Operator(1, AID(1, 1), AID(1, 2))
    @test AID(1, 1)*Operator(2, AID(1, 2)) == Operator(2, AID(1, 1), AID(1, 2)) == Operator(2, AID(1, 1))*AID(1, 2)
    @test Operator(2, AID(1, 1))//3 == Operator(2//3, AID(1, 1))

    @test +opts == opts
    @test -opts == Operators(-opt₁, -opt₂)
    @test opts*2 == 2*opts == Operators(2opt₁, 2opt₂)
    @test opts/2 == Operators(opt₁/2, opt₂/2)
    @test opts^2 == opts*opts == Operators(opt₁*opt₁, opt₁*opt₂, opt₂*opt₁, opt₂*opt₂)
    @test opts+1 == 1+opts == Operators(Operator(1), opt₁, opt₂)
    @test 1-opts == Operators(Operator(1), -opt₁, -opt₂)
    @test opts-1 == Operators(opt₁, opt₂, Operator(-1))
    @test opt₁+opt₂ == opts
    @test opt₁+opts == opts+opt₁ == Operators(opt₁*2, opt₂)
    @test opts+opts == Operators(opt₁*2, opt₂*2)
    @test opts-opt₁ == Operators(opt₂)
    @test opt₁-opts == Operators(-opt₂)
    @test opts-opts == zero(opts) == zero(typeof(opts))
    @test opts*opt == Operators(opt₁*opt, opt₂*opt)
    @test opt*opts == Operators(opt*opt₁, opt*opt₂)
end

@testset "LaTeX" begin
    @test latexname(AID) == :AID
    @test latexformat(AID) == LaTeX{(:nambu,), (:orbital,)}('c')

    latex = LaTeX{(:nambu,), (:orbital,)}('d')
    @test superscript(latex|>typeof) == (:nambu,)
    @test subscript(latex|>typeof) == (:orbital,)
    latexformat(AID, latex)

    aid = AID(1, 2)
    @test script(aid, latex, Val(:BD)) == 'd'
    @test script(aid, latex, Val(:SP)) == ("\\dagger",)
    @test script(aid, latex, Val(:SB)) == (1,)
    @test script(aid, Val(:spin)) == ""
    @test latexstring(aid) == "d^{\\dagger}_{1}"

    opt = Operator(1.0, AID(2, 2), AID(1, 1))
    @test latexstring(opt) == "d^{\\dagger}_{2}d^{}_{1}"
    io = IOBuffer()
    show(io, MIME"text/latex"(), opt)
    @test String(take!(io)) == "\$d^{\\dagger}_{2}d^{}_{1}\$"
    @test latexstring(Operator(1.0, aid, aid)) == "(d^{\\dagger}_{1})^2"

    latexformat(AID, LaTeX{(:nambu,), (:orbital,)}('c'))
    opts = Operators(
            Operator(1.0-1.0im, AID(2, 2), AID(1, 1)),
            Operator(-1.0, AID(1, 2), AID(1, 1))
            )
    candidate₁ = "(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1}"
    candidate₂ = "-c^{\\dagger}_{1}c^{}_{1}+(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}"
    @test latexstring(opts) ∈ (candidate₁, candidate₂)

    show(io, MIME"text/latex"(), [opts, opts])
    @test String(take!(io)) == "\\begin{bmatrix}(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1} & (1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1}\\end{bmatrix}"
    show(io, MIME"text/latex"(), [opts opts; opts opts])
    @test String(take!(io)) == "\\begin{pmatrix}(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1} & (1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1} \\\\ (1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1} & (1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1}\\end{pmatrix}"
end

struct DoubleCoeff <: LinearTransformation end
@inline Base.valtype(::Type{DoubleCoeff}, M::Type{<:Union{Operator, Operators}}) = M
@inline (double::DoubleCoeff)(m::Operator) = replace(m, value(m)*2)

@testset "LinearTransformation" begin
    m = Operator(1, AID(1, 1))
    s = Operators(m)

    i = LinearFunction(identity)
    @test i==deepcopy(i) && isequal(i, deepcopy(i))
    @test i(m)==m && i(s)==s

    double = DoubleCoeff()
    @test valtype(double, m) == valtype(typeof(double), typeof(m)) == typeof(m)
    @test valtype(double, s) == valtype(typeof(double), typeof(s)) == typeof(s)
    @test map!(double, s) == s == Operators(Operator(2, AID(1, 1)))
end

@testset "Permutation" begin
    id₁, id₂ = AID(1, 1), AID(1, 2)
    opt₀ = Operator(1.5)
    opt₁ = Operator(1.5, id₁, id₂)
    opt₂ = Operator(1.5, id₂, id₁)

    p = Permutation(Dict(id₁=>1, id₂=>2))
    @test p(opt₁) == opt₂-opt₀
    @test p(opt₁+opt₂) == 2*opt₂-opt₀

    p = Permutation(Dict(id₁=>2, id₂=>1))
    @test p(opt₂) == opt₁+opt₀
    @test p(opt₁+opt₂) == 2*opt₁+opt₀
end

@testset "Substitution" begin
    id₁, id₂ = AID(1, 1), AID(1, 2)
    opt = Operator(1.5, id₁, id₂)

    id₃, id₄ = AID(2, 1), AID(2, 2)
    ops₁ = Operator(2.0, id₃) + Operator(3.0, id₄) + 1
    ops₂ = Operator(2.0, id₃) - Operator(3.0, id₄)

    M = promote_type(eltype(ops₁), eltype(ops₂))
    table = Dict{eltype(opt), Operators{M, idtype(M)}}(id₁=>ops₁, id₂=>ops₂)
    substitution = TabledUnitSubstitution(table)
    @test substitution(opt) == 1.5*ops₁*ops₂
    @test substitution(Operators(opt)) == 1.5*ops₁*ops₂
end

@testset "RankFilter" begin
    op₀, op₁, op₂ = Operator(1), Operator(2, AID(1, 1)), Operator(3, AID(1, 1), AID(2, 1))
    ops = Operators(op₀, op₁, op₂)
    @test RankFilter(0)(ops) == Operators(op₀)
    @test RankFilter(1)(ops) == Operators(op₁)
    @test RankFilter(2)(ops) == Operators(op₂)
end
