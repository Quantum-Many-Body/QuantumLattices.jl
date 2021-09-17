using Test
using Printf: @printf
using QuantumLattices.Essentials.QuantumAlgebras
using QuantumLattices.Interfaces: id, value, rank, add!, sub!, mul!, div!
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.Combinatorics: Combinations
using QuantumLattices.Prerequisites.Traits: contentnames, parameternames, parametertype, isparameterbound

import QuantumLattices.Interfaces: ⊗, ⋅, permute

struct AID{O<:Real, S<:Real} <: SimpleID
    orbital::O
    nambu::S
end

@testset "ID" begin
    @test ID{SimpleID, 1}|>rank == 1
    @test ID|>rank == Any
    @test promote_type(ID{AID{Int, Int}, 1}, ID{AID{Int, Int}}) == ID{AID{Int, Int}}
    @test promote_type(ID{AID{Int, Int}}, ID{AID{Int, Int}, 1}) == ID{AID{Int, Int}}
    @test promote_type(ID{AID{Int, Int}}, ID{AID{Int, Float}}) == ID{AID{Int, <:Real}}
    @test promote_type(Tuple{}, ID{AID{Int, Int}}) == ID{AID{Int, Int}}
    @test promote_type(Tuple{}, ID{AID{Int, Int}, 2}) == ID{AID{Int, Int}}

    cid = AID(2, 1) * AID(1, Inf)
    @test cid == ID(AID(2, 1), AID(1, Inf))
    @test cid|>string == "ID(AID(2, 1), AID(1, Inf))"
    @test cid|>eltype == AID{Int, <:Real}
    @test cid|>rank == 2
    @test cid|>typeof|>rank == 2

    sid = AID(1, 1)
    @test sid * cid == ID(AID(1, 1), AID(2, 1), AID(1, Inf))
    @test cid * sid == ID(AID(2, 1), AID(1, Inf), AID(1, 1))
    @test cid * cid == ID(AID(2, 1), AID(1, Inf), AID(2, 1), AID(1, Inf))

    @test deepcopy(ID(AID(1, 1))) == ID(AID(1, 1))
    @test isless(SimpleID, ID(AID(1, 2)), ID(AID(1, Inf)))
    @test isless(SimpleID, ID(AID(2, 1)), ID(AID(1, 2), AID(1, Inf)))

    cid = AID(2, 1) * AID(3, 4.0)
    @test cid == ID(AID, (2, 3), (1, 4.0))
    @test cid|>propertynames == (:orbitals, :nambus)
    @test cid.orbitals == (2, 3)
    @test cid.nambus == (1, 4.0)
end

struct AOperator{V<:Number, I<:ID{SimpleID}} <: Element{V, I}
    value::V
    id::I
end
⊗(m₁::AOperator, m₂::AOperator) = m₁ * m₂
⋅(m₁::AOperator, m₂::AOperator) = m₁ * m₂
Base.show(io::IO, opt::AOperator) = @printf io "AOperator(value=%s, id=%s)" opt.value opt.id
Base.repr(opt::AOperator) = "AOperator($(opt.value), $(opt.id))"

function permute(id₁::AID, id₂::AID)
    @assert id₁ ≠ id₂ "permute error: permuted ids should not be equal to each other."
    if (id₁.nambu == 3-id₂.nambu) && (id₁.orbital == id₂.orbital)
        if id₁.nambu == 2
            return (AOperator(1, ID()), AOperator(1, ID(id₂, id₁)))
        else
            return (AOperator(-1, ID()), AOperator(1, ID(id₂, id₁)))
        end
    else
        return (AOperator(1, ID(id₂, id₁)),)
    end
end

@testset "Elements" begin
    @test contentnames(Element) == (:value, :id)
    @test parameternames(Element) == (:value, :id)
    @test isparameterbound(Element, :value, Any) == false
    @test isparameterbound(Element, :id, ID) == true
    @test isparameterbound(Element, :id, ID{AID{Int, Int}, 2}) == false
    @test valtype(Element) == parametertype(Element, :value) == parametertype(Element, 1) == Any
    @test valtype(Element{Int}) == parametertype(Element{Int}, :value) == parametertype(Element{Int}, 1) == Int
    @test idtype(Element{<:Number}) == parametertype(Element{<:Number}, :id) == parametertype(Element{<:Number}, 2) == ID{SimpleID}
    @test idtype(Element{<:Number, ID{AID}}) == parametertype(Element{<:Number, ID{AID}}, :id) == parametertype(Element{<:Number, ID{AID}},2) == ID{AID}
    @test promote_type(AOperator{Int}, AOperator) == AOperator
    @test promote_type(AOperator, AOperator{Int}) == AOperator
    @test promote_type(AOperator{Int}, AOperator{Float}) == AOperator{Float}
    @test promote_type(AOperator{Int, ID{AID{Int, Int}, 2}}, AOperator{Float, ID{AID{Int, Int}, 2}}) == AOperator{Float, ID{AID{Int, Int}, 2}}
    @test promote_type(AOperator{Int, ID{AID}}, AOperator{Float, ID{AID}}) == AOperator{Float, <:ID{AID}}

    opt = AOperator(1.0, ID(AID(1, 1)))
    @test value(opt) == 1.0
    @test id(opt) == ID(AID(1, 1))
    @test opt|>deepcopy == opt
    @test isequal(opt|>deepcopy, opt)
    @test isapprox(opt, replace(opt, value=opt.value+10^-6); atol=10^-5)
    @test opt|>valtype == opt|>typeof|>valtype == parametertype(opt|>typeof, :value) == Float
    @test opt|>idtype == opt|>typeof|>idtype == parametertype(opt|>typeof, :id) == ID{AID{Int, Int}, 1}
    @test opt|>rank == opt|>typeof|>rank == 1
    @test opt[1] == AOperator(1, ID(AID(1, 1)))
    @test length(opt) == 1
    @test firstindex(opt) == 1
    @test lastindex(opt) ==1
    @test opt|>typeof|>one == opt|>one == AOperator(1.0, ID())
    @test +opt == opt
    @test opt+zero(Elements) == zero(Elements)+opt == opt-zero(Elements)
    @test -opt == AOperator(-1.0, ID(AID(1, 1)))
    @test zero(Elements)-opt == -opt
    @test opt*2 == 2*opt == opt*AOperator(2, ID()) == AOperator(2, ID())*opt == AOperator(2.0, ID(AID(1, 1)))
    @test opt*zero(Elements) == zero(Elements)*opt == zero(Elements) == nothing
    @test opt/2 == opt/AOperator(2, ID()) == AOperator(0.5, ID(AID(1, 1)))
    @test opt^2 == opt*opt

    @test AOperator(1, ID())+1 == 1+AOperator(1, ID()) == AOperator(1, ID())+AOperator(1, ID()) == AOperator(2, ID())
    @test opt+1 == 1+opt == opt+AOperator(1, ID()) == AOperator(1, ID())+opt
    @test AOperator(1, ID())-2 == 1-AOperator(2, ID()) == AOperator(1, ID())-AOperator(2, ID()) == AOperator(-1, ID())
    @test AOperator(2, ID())*AOperator(3, ID()) == AOperator(6, ID())

    opt₁ = AOperator(1.0, ID(AID(1, 1)))
    opt₂ = AOperator(2.0, ID(AID(1, 2)))
    opts = Elements(opt₁.id=>opt₁, opt₂.id=>opt₂)
    @test string(opts) == "Elements with 2 entries:\n  AOperator(value=2.0, id=ID(AID(1, 2)))\n  AOperator(value=1.0, id=ID(AID(1, 1)))\n"
    @test repr(opts) == "Elements with 2 entries:\n  AOperator(2.0, ID(AID(1, 2)))\n  AOperator(1.0, ID(AID(1, 1)))"
    @test opts|>deepcopy == opts
    @test opts|>zero == opts|>typeof|>zero == nothing
    @test isapprox(opts, deepcopy(opts)) && isapprox(opts, opts+10^-6; atol=10^-5) && isapprox(opts+10^-6, opts; atol=10^-5)
    @test !isapprox(opts, opts+10^-6) && !isapprox(opts+10^-6, opts)
    @test add!(deepcopy(opts)) == opts
    @test add!(deepcopy(opts), zero(Elements)) == opts
    @test sub!(deepcopy(opts)) == opts
    @test sub!(deepcopy(opts), zero(Elements)) == opts
    @test mul!(deepcopy(opts), 2.0) == mul!(deepcopy(opts), AOperator(2, ID())) == opts*2 == opts*AOperator(2, ID())
    @test div!(deepcopy(opts), 2.0) == div!(deepcopy(opts), AOperator(2, ID())) == opts/2 == opts/AOperator(2, ID())
    @test +opts == opts+zero(Elements) == zero(Elements)+opts == opts == opts-zero(Elements)
    @test opt₁+opt₂ == opts
    @test opt₁+opts == opts+opt₁ == Elements(opt₁*2, opt₂)
    @test opts+opts == Elements(opt₁*2, opt₂*2)
    @test -opts == Elements(-opt₁, -opt₂) == zero(Elements)-opts
    @test 1.0-opts == AOperator(1.0, ID())-opts == -opts+1.0
    @test opts-1.0 == opts-AOperator(1.0, ID())
    @test opt-1.0 == opt-AOperator(1.0, ID()) == opt+AOperator(-1.0, ID())
    @test 1.0-opt == AOperator(1.0, ID())-opt == 1.0+(-opt)
    @test zero(Elements) == opt₁-opt₁ == opts-opts
    @test opts-opts == zero(Elements)
    @test isequal(zero(Elements), opts-opts)
    @test isequal(opts-opts, zero(Elements))
    @test opts-opt₁ == Elements(opt₂)
    @test opt₁-opts == Elements(-opt₂)
    @test 2*opts == opts*2 == Elements(opt₁*2, opt₂*2) == AOperator(2, ID())*opts
    @test opts/2 == Elements(opt₁/2, opt₂/2)
    @test opts^2 == opts*opts
    @test opts*zero(Elements) == zero(Elements)*opts == nothing

    temp = Elements{ID{AID{Int, Int}}, AOperator{Float}}(opts)
    @test add!(deepcopy(temp), 1) == add!(deepcopy(temp), AOperator(1, ID())) == opts+1 == 1+opts == opts+AOperator(1, ID()) == AOperator(1, ID())+opts
    @test sub!(deepcopy(temp), 1) == sub!(deepcopy(temp), AOperator(1, ID())) == opts-1 == opts-AOperator(1, ID())

    opt₃ = AOperator(3, ID(AID(1, 3)))
    @test opt₁*opt₂ == opt₁⊗opt₂ == opt₁⋅opt₂ == AOperator(2.0, ID(AID(1, 1), AID(1, 2)))
    @test opts*opt₃ == opts⊗opt₃ == opts⋅opt₃ == Elements(AOperator(3.0, ID(AID(1, 1), AID(1, 3))), AOperator(6.0, ID(AID(1, 2), AID(1, 3))))
    @test opt₃*opts == opt₃⊗opts == opt₃⋅opts == Elements(AOperator(3.0, ID(AID(1, 3), AID(1, 1))), AOperator(6.0, ID(AID(1, 3), AID(1, 2))))
    @test opts*Elements(opt₃) == opts⊗Elements(opt₃) == opts⋅Elements(opt₃) == Elements(AOperator(3.0, ID(AID(1, 1), AID(1, 3))), AOperator(6.0, ID(AID(1, 2), AID(1, 3))))

    table = Dict(AID(1, 1)=>1, AID(1, 2)=>2)
    opt = AOperator(2.0, ID(AID(1, 1), AID(1, 2)))
    @test sequence(opt, table) == (1, 2)
    @test split(opt) == (2.0, AOperator(1.0, ID(AID(1, 1))), AOperator(1.0, ID(AID(1, 2))))

    opt₁ = AOperator(1, ID(AID(1, 1)))
    opt₂ = AOperator(2, ID(AID(1, 2)))
    opts = Elements(opt₁.id=>opt₁, opt₂.id=>opt₂)
    opt = AOperator(3, ID())
    @test opt₁//3 == opt₁//opt == AOperator(1//3, ID(AID(1, 1)))
    @test opt₂//3 == opt₂//opt == AOperator(2//3, ID(AID(1, 2)))
    @test opts//3 == opts//opt == Elements(opt₁.id=>AOperator(1//3, ID(AID(1, 1))), opt₂.id=>AOperator(2//3, ID(AID(1, 2))))
end

@testset "replace" begin
    id₁, id₂ = AID(1, 1), AID(1, 2)
    opt = AOperator(1.5, ID(id₁, id₂))
    id₃, id₄ = AID(2, 1), AID(2, 2)
    opt₁ = AOperator(2.0, ID(id₃))+AOperator(3.0, ID(id₄))
    opt₂ = AOperator(2.0, ID(id₃))-AOperator(3.0, ID(id₄))
    @test replace(opt, id₁=>opt₁, id₂=>opt₂) == 1.5*opt₁*opt₂
    @test replace(Elements(opt), id₂=>opt₂, id₁=>opt₁) == 1.5*opt₁*opt₂
end

@testset "permute" begin
    id₁, id₂ = AID(1, 1), AID(1, 2)
    opt₁ = AOperator(1.5, ID(id₁, id₂))
    opt₂ = AOperator(1.5, ID(id₂, id₁))
    opt₀ = AOperator(1.5, ID())
    @test permute(opt₁, Dict(id₁=>1, id₂=>2)) == opt₂-opt₀
    @test permute(opt₁+opt₂, Dict(id₁=>1, id₂=>2)) == 2*opt₂-opt₀
    @test permute(opt₂, Dict(id₁=>2, id₂=>1)) == opt₁+opt₀
    @test permute(opt₁+opt₂, Dict(id₁=>2, id₂=>1)) == 2*opt₁+opt₀
end

@testset "Transformation" begin
    m = AOperator(1.0, ID(AID(1, 1)))
    s = Elements(m)
    i = Identity()
    @test valtype(i, m) == valtype(typeof(i), typeof(m)) == typeof(m)
    @test valtype(i, s) == valtype(typeof(i), typeof(s)) == typeof(s)
    @test i(m)==m && i(s)==s
end
