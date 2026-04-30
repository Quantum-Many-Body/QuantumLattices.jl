using Latexify: latexify
using LinearAlgebra: dot
using Printf: @sprintf
using QuantumLattices: ZeroAtLeast, add!, div!, id, mul!, rank, sub!, value, ⊗
using QuantumLattices.QuantumOperators
using QuantumLattices.Toolkit: Float, contentnames, isparameterbound, parameternames, parametertype, subscript, superscript

import QuantumLattices: permute
import QuantumLattices.QuantumOperators: script

# Bose
struct Bose{O<:Real, S<:Real} <: OperatorIndex
    orbital::O
    nambu::S
end
@inline Base.adjoint(id::Bose) = Bose(id.orbital, 3-id.nambu)
@inline script(id::Bose, ::Val{:orbital}; kwargs...) = id.orbital
@inline script(id::Bose, ::Val{:nambu}; kwargs...) = id.nambu==2 ? "\\dagger" : ""
latexformat(Bose, LaTeX{(:nambu,), (:orbital,)}('c'))
function permute(u₁::Bose, u₂::Bose)
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
    @test promote_type(ZeroAtLeast{Bose{Int, Int}, 1}, ZeroAtLeast{Bose{Int, Int}}) == ZeroAtLeast{Bose{Int, Int}}
    @test promote_type(ZeroAtLeast{Bose{Int, Int}}, ZeroAtLeast{Bose{Int, Int}, 1}) == ZeroAtLeast{Bose{Int, Int}}
    @test promote_type(ZeroAtLeast{Bose{Int, Int}}, ZeroAtLeast{Bose{Int, Float}}) == ZeroAtLeast{Bose{Int, <:Real}}
    @test promote_type(Tuple{}, ZeroAtLeast{Bose{Int, Int}}) == ZeroAtLeast{Bose{Int, Int}}
    @test promote_type(Tuple{}, ZeroAtLeast{Bose{Int, Int}, 2}) == ZeroAtLeast{Bose{Int, Int}}

    sid = Bose(1, 1)
    @test string(sid) == "Bose(1, 1)"
    @test !iszero(sid)
    @test sid'==Bose(1, 2)
    @test hash(sid) == hash((1, 1))
    @test OperatorIndex[sid] == "Bose"

    cid = (Bose(2, 1), Bose(Inf, 1))
    @test cid == (Bose(2, 1) ⊗ Bose(Inf, 1))
    @test ⊗(sid, cid) == ⊗(Bose(1, 1), Bose(2, 1), Bose(Inf, 1))
    @test ⊗(cid, sid) == ⊗(Bose(2, 1), Bose(Inf, 1), Bose(1, 1))
    @test ⊗(cid, cid) == ⊗(Bose(2, 1), Bose(Inf, 1), Bose(2, 1), Bose(Inf, 1))
    @test cid == ZeroAtLeast(Bose, (2, Inf), (1, 1))
    @test cid|>propertynames == (:orbitals, :nambus)
    @test cid.orbitals == (2, Inf)
    @test cid.nambus == (1, 1)
    @test cid|>eltype == Bose{<:Real, Int}
    @test cid|>rank == cid|>typeof|>rank == 2
    @test cid'==(Bose(Inf, 2), Bose(2, 2))
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
    @test idtype(OperatorProd{<:Number, ZeroAtLeast{Bose}}) == parametertype(OperatorProd{<:Number, ZeroAtLeast{Bose}}, :id) == parametertype(OperatorProd{<:Number, ZeroAtLeast{Bose}}, 2) == ZeroAtLeast{Bose}

    @test promote_type(Operator{Int}, Operator) == Operator
    @test promote_type(Operator, Operator{Int}) == Operator
    @test promote_type(Operator{Int}, Operator{Float}) == Operator{Float}
    @test promote_type(Operator{Int, ZeroAtLeast{Bose{Int, Int}, 2}}, Operator{Float, ZeroAtLeast{Bose{Int, Int}, 2}}) == Operator{Float, ZeroAtLeast{Bose{Int, Int}, 2}}
    @test promote_type(Operator{Int, ZeroAtLeast{Bose}}, Operator{Float, ZeroAtLeast{Bose}}) == Operator{Float, <:ZeroAtLeast{Bose}}
    @test promote_type(Operator{Int}, Float) == Operator{Float}

    opt = Operator(2.0, Bose(1, 1))
    @test deepcopy(opt)==opt && isequal(deepcopy(opt), opt)
    @test eltype(opt) == eltype(typeof(opt)) == Bose{Int, Int}
    @test collect(opt) == [Bose(1, 1)]
    @test opt|>rank == opt|>typeof|>rank == 1
    @test opt|>valtype == opt|>typeof|>valtype == parametertype(opt|>typeof, :value) == Float
    @test opt|>idtype == opt|>typeof|>idtype == parametertype(opt|>typeof, :id) == ZeroAtLeast{Bose{Int, Int}, 1}
    @test opt|>scalartype == opt|>typeof|>scalartype == Float
    @test opt|>isequivalenttoscalar == opt|>typeof|>isequivalenttoscalar == false
    @test value(opt) == 2.0
    @test id(opt) == ⊗(Bose(1, 1))
    @test replace(opt, 3) == Operator(3, Bose(1, 1))
    @test !iszero(opt) && iszero(replace(opt, 0))
    @test isapprox(opt, replace(opt, opt.value+10^-6); atol=10^-5)
    @test length(opt) == 1 && firstindex(opt) == 1 && lastindex(opt) ==1
    @test opt[1]==opt[begin]==opt[end]==Bose(1, 1) && opt[1:1]==Operator(1.0, Bose(1, 1))
    @test split(opt) == (2.0, Bose(1, 1))
    @test opt|>typeof|>one == opt|>one == Operator(1.0)
    @test convert(Operator{Float, Tuple{Bose{Int, Int}}}, Operator(2, Bose(1, 1))) == Operator(2.0, Bose(1, 1))
    @test Operator(2.0) == Operator(2.0)
    @test string(opt) == "Operator(2.0, Bose(1, 1))"
    @test opt' == Operator(2.0, Bose(1, 2))
    @test ishermitian(opt) == false
    @test ishermitian(Operator(2.0, Bose(1, 1), Bose(1, 2)))
    @test sequence(opt, Dict(Bose(1, 1)=>1, Bose(1, 2)=>2)) == (1,)
    @test convert(typeof(opt), Bose(1, 1)) == Operator(1.0, Bose(1, 1))
    @test convert(Operator{Float, <:ZeroAtLeast{Bose}}, Bose(1, 1)) == Operator(1.0, Bose(1, 1))
    @test convert(Operator{Float, Tuple{}}, 2) == Operator(2.0)
    @test convert(Operator{Float, <:ZeroAtLeast{Bose}}, 2) == Operator(2.0)

    opt₁ = Operator(1.0, Bose(1, 1))
    opt₂ = Operator(2.0, Bose(1, 2))
    opts = Operators(opt₁, opt₂)
    @test opts == Operators{eltype(opts)}(opt₁, opt₂) == OperatorSum(opt₁, opt₂) == OperatorSum((opt₁, opt₂)) == OperatorSum{eltype(opts)}(opt₁, opt₂) == OperatorSum{eltype(opts)}((opt₁, opt₂))
    @test eltype(opts) == eltype(typeof(opts)) == Operator{Float, ZeroAtLeast{Bose{Int, Int}, 1}}
    @test scalartype(opts) == scalartype(typeof(opts)) == Float
    @test isequivalenttoscalar(opts) == isequivalenttoscalar(typeof(opts)) == false
    @test collect(opts) == collect(values(opts.contents))
    @test length(opts) == 2
    @test summary(opts) == "Operators"
    @test string(opts) == @sprintf "Operators with 2 Operator\n  %s\n" join(opts, "\n  ")
    @test haskey(opts, id(opt₁)) && haskey(opts, id(opt₂)) && !haskey(opts, ⊗(Bose(3, 1)))
    @test opts[1]==opts[begin]==opt₁ && opts[2]==opts[end]==opt₂
    @test opts[1:2] == opts[:] == opts
    @test empty(opts) == empty!(deepcopy(opts)) == zero(opts)
    @test !iszero(opts) && iszero(zero(opts))
    @test !invoke(iszero, Tuple{OperatorSet}, opts) && invoke(iszero, Tuple{OperatorSet}, zero(opts))
    optp₁ = promote_type(typeof(opts), OperatorSum{Operator{Complex{Int}, NTuple{2, Bose{Float, Float}}}, NTuple{2, Bose{Float, Float}}})
    optp₂ = OperatorSum{Operator{Complex{Float}, <:Tuple{Bose, Vararg{Bose{Float, Float}}}}, Tuple{Bose, Vararg{Bose{Float, Float}}}}
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

    @test operatortype(Bose(1, 1)) == operatortype(Bose{Int, Int}) == Operator{Int, Tuple{Bose{Int, Int}}}
    @test operatortype(Bose) == Operator{Int, <:Tuple{Bose}}
    @test operatortype(opt) == operatortype(typeof(opt)) == Operator{Float, Tuple{Bose{Int, Int}}}
    @test operatortype(opts) == operatortype(typeof(opts)) == eltype(opts)

    @test zero(Bose(1, 1)) == zero(Operator(1, Bose(1, 1))) == zero(Operator(1, Bose(1, 1))+Operator(1, Bose(2, 1)))
    @test conj(Bose(1, 1)) == Bose(1, 1)
    @test conj(Operator(2im, Bose(2, 1), Bose(1, 1))) == Operator(-2im, Bose(2, 1), Bose(1, 1))
    @test conj(Operator(2im, Bose(2, 1), Bose(1, 1))+Operator(2, Bose(1, 1), Bose(1, 1))) == Operator(-2im, Bose(2, 1), Bose(1, 1))+Operator(2, Bose(1, 1), Bose(1, 1))
    @test dot(Operator(2im, Bose(1, 1)), Operator(2im, Bose(1, 1), Bose(2, 1))) == Operator(4, Bose(1, 1), Bose(1, 1), Bose(2, 1))
    @test dot(Operator(2im, Bose(2, 1), Bose(1, 1)), 2) == Operator(-4im, Bose(2, 1), Bose(1, 1))
    @test dot(2im, Operator(2im, Bose(2, 1), Bose(1, 1))) == Operator(4, Bose(2, 1), Bose(1, 1))

    @test +opt == opt
    @test -opt == Operator(-2.0, Bose(1, 1))
    @test opt*2 == 2*opt == Operator(4.0, Bose(1, 1))
    @test opt/2 == Operator(1.0, Bose(1, 1))
    @test opt^2 == opt*opt == Operator(4.0, Bose(1, 1), Bose(1, 1))
    @test opt+1 == 1+opt == Operators(Operator(1), opt)
    @test 1-opt == Operators(Operator(1), -opt)
    @test opt-1 == Operators(opt, Operator(-1))
    @test 2*Bose(1, 1) == Bose(1, 1)*2 == Operator(2, Bose(1, 1))
    @test Bose(1, 1)*Bose(1, 2) == Operator(1, Bose(1, 1), Bose(1, 2))
    @test Bose(1, 1)*Operator(2, Bose(1, 2)) == Operator(2, Bose(1, 1), Bose(1, 2)) == Operator(2, Bose(1, 1))*Bose(1, 2)
    @test Operator(2, Bose(1, 1))//3 == Operator(2//3, Bose(1, 1))

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
    @test latexname(Bose) == :Bose
    @test latexformat(Bose) == LaTeX{(:nambu,), (:orbital,)}('c')

    latex = LaTeX{(:nambu,), (:orbital,)}('d')
    @test superscript(latex|>typeof) == (:nambu,)
    @test subscript(latex|>typeof) == (:orbital,)
    latexformat(Bose, latex)

    aid = Bose(1, 2)
    @test script(aid, latex, Val(:BD)) == 'd'
    @test script(aid, latex, Val(:SP)) == ("\\dagger",)
    @test script(aid, latex, Val(:SB)) == (1,)
    @test script(aid, Val(:spin)) == ""
    @test String(latexify(aid; env=:raw)) == "d^{\\dagger}_{1}"

    opt = Operator(1.0, Bose(2, 2), Bose(1, 1))
    @test String(latexify(opt; env=:raw)) == "d^{\\dagger}_{2}d^{}_{1}"
    io = IOBuffer()
    show(io, MIME"text/latex"(), opt)
    @test String(take!(io)) == "\$d^{\\dagger}_{2}d^{}_{1}\$"
    @test String(latexify(Operator(1.0, aid, aid); env=:raw)) == "(d^{\\dagger}_{1})^2"

    latexformat(Bose, LaTeX{(:nambu,), (:orbital,)}('c'))
    opts = Operators(
            Operator(1.0-1.0im, Bose(2, 2), Bose(1, 1)),
            Operator(-1.0, Bose(1, 2), Bose(1, 1))
            )
    str = "\\left(1.0-1.0\\mathit{i}\\right)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1}"
    @test String(latexify(opts; env=:raw)) == str

    show(io, MIME"text/latex"(), [opts, opts])
    @test String(take!(io)) == "\\begin{equation}\n\\left[\n\\begin{array}{c}\n$str \\\\\n$str \\\\\n\\end{array}\n\\right]\n\\end{equation}\n"
    show(io, MIME"text/latex"(), [opts opts; opts opts])
    @test String(take!(io)) == "\\begin{equation}\n\\left[\n\\begin{array}{cc}\n$str & $str \\\\\n$str & $str \\\\\n\\end{array}\n\\right]\n\\end{equation}\n"
end

struct DoubleCoeff <: LinearTransformation end
@inline Base.valtype(::Type{DoubleCoeff}, M::Type{<:Union{Operator, Operators}}) = M
@inline (double::DoubleCoeff)(m::Operator) = replace(m, value(m)*2)

@testset "LinearTransformation" begin
    m = Operator(1, Bose(1, 1))
    s = Operators(m)

    double = DoubleCoeff()
    @test valtype(double, m) == valtype(typeof(double), typeof(m)) == typeof(m)
    @test valtype(double, s) == valtype(typeof(double), typeof(s)) == typeof(s)
    @test map!(double, s) == s == Operators(Operator(2, Bose(1, 1)))

    i = LinearFunction(identity)
    @test i==deepcopy(i) && isequal(i, deepcopy(i))
    @test i(m)==m && i(s)==s
    @test valtype(LinearFunction{typeof(m), typeof(identity)}, typeof(m)) == typeof(m)
    @test valtype(LinearFunction{typeof(m), typeof(identity)}, typeof(s)) == typeof(s)

    op₀, op₁, op₂ = Operator(1), Operator(2, Bose(1, 1)), Operator(3, Bose(1, 1), Bose(2, 1))
    ops = Operators(op₀, op₁, op₂)
    @test LinearFunction{Operator{Int, Tuple{}}}(op->rank(op)==0 ? op : 0)(ops) == Operators(op₀)
    @test LinearFunction{Operator{Int, Tuple{Bose{Int, Int}}}}(op->rank(op)==1 ? op : 0)(ops) == Operators(op₁)
    @test LinearFunction{Operator{Int, NTuple{2, Bose{Int, Int}}}}(op->rank(op)==2 ? op : 0)(ops) == Operators(op₂)
end

@testset "Permutation" begin
    id₁, id₂ = Bose(1, 1), Bose(1, 2)
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
    id₁, id₂ = Bose(1, 1), Bose(1, 2)
    opt = Operator(1.5, id₁, id₂)

    id₃, id₄ = Bose(2, 1), Bose(2, 2)
    ops₁ = Operator(2.0, id₃) + Operator(3.0, id₄) + 1
    ops₂ = Operator(2.0, id₃) - Operator(3.0, id₄)

    M = promote_type(eltype(ops₁), eltype(ops₂))
    table = Dict{eltype(opt), Operators{M, idtype(M)}}(id₁=>ops₁, id₂=>ops₂)
    substitution = TabledUnitSubstitution(table)
    @test substitution(opt) == 1.5*ops₁*ops₂
    @test substitution(Operators(opt)) == 1.5*ops₁*ops₂
end
