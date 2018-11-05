using Test
import MacroTools: rmlines
using Hamiltonian.Utilities.Factory

rmlines(exprs::Vector)=filter(expr->!isa(expr,LineNumberNode),exprs)

@testset "Block" begin
    block=Block(:(x=1;y=2))
    @test block.body|>rmlines==Any[:(x=1),:(y=2)]
    @test block()|>rmlines==Expr(:block,:(x=1),:(y=2))

    block=@block a=1
    @push! block b=2 c=3 d=4
    @test block.body==[:(a=1),:(b=2),:(c=3),:(d=4)]
end

@testset "Argument" begin
    argument=Argument(:(phone::Vector{Int}=[9,1,1]))
    @test argument.name==:phone
    @test argument.type==:(Vector{Int})
    @test argument.slurp==false
    @test argument.default==:([9,1,1])
    @test argument()==Expr(:kw,:(phone::Vector{Int}),:[9,1,1])

    argument=Argument(:(properties::Int...))
    @test argument.name==:properties
    @test argument.type==:Int
    @test argument.slurp==true
    @test argument.default==nothing
    @test argument()==:(properties::Int...)

    argument1=@argument address::String
    argument2=Argument(:(address::String))
    @test argument1==argument2
end

@testset "FunctionFactory" begin
    ff=FunctionFactory(:(fx(x::Int,y::Int;choice::Function=sum)=choice(x,y)))
    @test ff.name==:fx
    @test ff.args==Union{Symbol,Expr}[:(x::Int),:(y::Int)]
    @test ff.kwargs==Union{Symbol,Expr}[Expr(:kw,:(choice::Function),:sum)]
    @test ff.rtype==:Any
    @test ff.params==Union{Symbol,Expr}[]
    @test ff.body==:(begin choice(x,y) end)

    ff=FunctionFactory(:(
        function fx(x::T,y::T;choice::Function=sum) where T
            choice(x,y)
        end
    ))
    @test ff.name==:fx
    @test ff.args==Union{Symbol,Expr}[:(x::T),:(y::T)]
    @test ff.kwargs==Union{Symbol,Expr}[Expr(:kw,:(choice::Function),:sum)]
    @test ff.rtype==:Any
    @test ff.params==Union{Symbol,Expr}[:T]
    @test ff.body|>rmlines==:(begin choice(x,y) end)

    ff=@functionfactory (fx()::T) where T<:Number=nothing
    @addargs! ff x::T y::T z::T
    @addkwargs! ff sign::Int=1 choice::Function=sum
    @extendbody! ff result=choice(x,y,z) result*=sign
    @test ff.name==:fx
    @test ff.args==Union{Symbol,Expr}[:(x::T),:(y::T),:(z::T)]
    @test ff.kwargs==Union{Symbol,Expr}[Expr(:kw,:(sign::Int),:1),Expr(:kw,:(choice::Function),:sum)]
    @test ff.rtype==:T
    @test ff.params==Union{Symbol,Expr}[:(T<:Number)]
    @test ff.body|>rmlines==Expr(:block,:nothing,:(result=choice(x,y,z)),:(result*=sign))
end

@testset "Field" begin
    field=Field(:xs)
    @test field.name==:xs
    @test field.type==:Any
    @test field()==:(xs::Any)

    field=Field(:(xs::Vector{Int}))
    @test field.name==:xs
    @test field.type==:(Vector{Int})
    @test field()==:(xs::Vector{Int})

    field1=@field xs::Vector{Int}
    field2=Field(:(xs::Vector{Int}))
    @test field1==field2
end

@testset "TypeFactory" begin
    tf=TypeFactory(:(
        struct Child{T} <: Parent{T} where T<:Real
            address::String
            phone::Vector{Int}
            properties::Vector{T}
        end
    ))
    @test tf.name==:Child
    @test tf.mutable==false
    @test tf.params==Union{Symbol,Expr}[:T]
    @test tf.supertype==:(Parent{T} where T <: Real)
    @test tf.fields==Tuple{Symbol,Union{Symbol,Expr}}[(:address,:String),(:phone,:(Vector{Int})),(:properties,:(Vector{T}))]
    @test tf.constructors==Expr[]

    tf=@typefactory struct Child{T} <: Parent{T} where T<:Real end
    @addfields! tf address::String phone::Vector{Int} properties::Vector{T}
    @addconstructors! tf Child(address::String,phone::Vector{Int},properties::T)=new(address,phone,properties)
    @test tf.name==:Child
    @test tf.mutable==false
    @test tf.params==Union{Symbol,Expr}[:T]
    @test tf.supertype==:(Parent{T} where T <: Real)
    @test tf.fields==Tuple{Symbol,Union{Symbol,Expr}}[(:address,:String),(:phone,:(Vector{Int})),(:properties,:(Vector{T}))]
end
