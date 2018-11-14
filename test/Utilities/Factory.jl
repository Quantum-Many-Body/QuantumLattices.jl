using Test
using Hamiltonian.Utilities.Factory

@testset "Inference" begin
    unescaped=(:Vector,:Tuple,:Int,:String)
    @test (@inference Vector{Tuple{Int,String}} unescaped)()==:(Vector{Tuple{Int,String}})
    unescaped=(:Type,:Real)
    @test (@inference Type{<:Real} unescaped)()==:(Type{<:Real})
end

@testset "Argument" begin
    unescaped=(:Vector,:Int)
    argument=Argument(:(phone::Vector{Int}=[9,1,1]),unescaped=unescaped)
    @test argument.name==:phone
    @test argument.type==Inference(:(Vector{Int}),unescaped=unescaped)
    @test argument.slurp==false
    @test argument.default==:([9,1,1])
    @test argument()==Expr(:kw,:(phone::Vector{Int}),:[9,1,1])

    unescaped=(:Int,)
    argument=Argument(:(properties::Int...),unescaped=unescaped)
    @test argument.name==:properties
    @test argument.type==Inference(:Int,unescaped=unescaped)
    @test argument.slurp==true
    @test argument.default==nothing
    @test argument()==:(properties::Int...)

    unescaped=(:String,)
    argument1=@argument address::String unescaped
    argument2=Argument(:(address::String),unescaped=unescaped)
    @test argument1==argument2
end

@testset "Parameter" begin
    unescaped=()
    parameter=Parameter(:T,unescaped=unescaped)
    @test parameter.name==:T
    @test parameter.type==nothing
    @test parameter()==:T

    unescaped=(:Real,)
    parameter=Parameter(:(T<:Real),unescaped=unescaped)
    @test parameter.name==:T
    @test parameter.type==Inference(:Real,unescaped=unescaped)
    @test parameter()==:(T<:Real)

    unescaped=(:Real,)
    parameter=Parameter(:(<:Real),unescaped=unescaped)
    @test parameter.name==nothing
    @test parameter.type==Inference(:Real,unescaped=unescaped)
    @test parameter()==:(<:Real)

    unescaped=(:Int,)
    parameter=@parameter ::Int unescaped=unescaped
    @test parameter.name==nothing
    @test parameter.type==Inference(:Int,unescaped=unescaped)
    @test parameter()==:(<:Int)
end

@testset "Field" begin
    unescaped=(:Any,)
    field=Field(:xs,unescaped=unescaped)
    @test field.name==:xs
    @test field.type==Inference(:Any,unescaped=unescaped)
    @test field()==:(xs::Any)

    unescaped=(:Vector,:Int)
    field=Field(:(xs::Vector{Int}),unescaped=unescaped)
    @test field.name==:xs
    @test field.type==Inference(:(Vector{Int}),unescaped=unescaped)
    @test field()==:(xs::Vector{Int})

    unescaped=(:Vector,:Int)
    field1=@field xs::Vector{Int} unescaped
    field2=Field(:(xs::Vector{Int}),unescaped=unescaped)
    @test field1==field2
end

@testset "Block" begin
    block=@rmlines Block(:(x=1;y=2))
    @test block.body==Any[:(x=1),:(y=2)]
    @test block()==Expr(:block,:(x=1),:(y=2))

    block=@block a=1
    @push! block b=2 c=3 d=4
    @test block.body==[:(a=1),:(b=2),:(c=3),:(d=4)]
end

@testset "FunctionFactory" begin
    unescaped=()
    ff=FunctionFactory(:(fx(x::Int,y::Int;choice::Function=sum)=choice(x,y)),unescaped=unescaped)
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::Int),:(y::Int)],unescaped=unescaped)
    @test ff.kwargs==Argument.([:(choice::Function=sum)],unescaped=unescaped)
    @test ff.rtype==Inference(:Any,unescaped=unescaped)
    @test ff.params==Parameter[]
    @test ff.body==Block(:(begin choice(x,y) end))

    unescaped=(:T,)
    ff=FunctionFactory(:(
        function fx(x::T,y::T;choice::Function=sum) where T
            choice(x,y)
        end
    ),unescaped=unescaped)
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::T),:(y::T)],unescaped=unescaped)
    @test ff.kwargs==Argument.([:(choice::Function=sum)],unescaped=unescaped)
    @test ff.rtype==Inference(:Any,unescaped=unescaped)
    @test ff.params==Parameter.([:T],unescaped=unescaped)
    @test ff.body|>rmlines==Block(:(begin choice(x,y) end))

    ff=@functionfactory (fx()::T)=nothing
    @addargs! ff x::T y::T z::T
    @addkwargs! ff sign::Int=1 choice::Function=sum
    @addparams! ff T<:Number
    @extendbody! ff result=choice(x,y,z) result*=sign
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::T),:(y::T),:(z::T)])
    @test ff.kwargs==Argument.([:(sign::Int=1),:(choice::Function=sum)])
    @test ff.rtype==Inference(:T)
    @test ff.params==Parameter.([:(T<:Number)])
    @test ff.body|>rmlines==Block(:nothing,:(result=choice(x,y,z)),:(result*=sign))
end

@testset "TypeFactory" begin
    unescaped=(:T,)
    tf=TypeFactory(:(
        struct Child{T<:Real} <: Parent{T}
            address::String
            phone::Vector{Int}
            properties::Vector{T}
        end
    ),unescaped=unescaped)
    @test tf.name==:Child
    @test tf.mutable==false
    @test tf.params==Parameter.([:(T<:Real)],unescaped=unescaped)
    @test tf.supertype==Inference(:(Parent{T}),unescaped=unescaped)
    @test tf.fields==Field.([:(address::String),:(phone::Vector{Int}),:(properties::Vector{T})],unescaped=unescaped)
    @test tf.constructors==FunctionFactory[]

    tf=@typefactory struct Child <: Parent{T} end
    @addparams! tf T<:Real
    @addfields! tf address::String phone::Vector{Int} properties::Vector{T}
    @addconstructors! tf Child(address::String,phone::Vector{Int},properties::Vector{T}) where T=new(address,phone,properties)
    @test tf.name==:Child
    @test tf.mutable==false
    @test tf.params==Parameter.([:(T<:Real)])
    @test tf.supertype==Inference(:(Parent{T}))
    @test tf.fields==Field.([:(address::String),:(phone::Vector{Int}),:(properties::Vector{T})])
end

eval((@typefactory struct Child name::String; seniority::Int; Child(name::String,seniority::Int=1)=new(name,seniority) end (:String,:Int,:Any))())

@testset "call" begin
    c=Child("Tuanzi")
    @test c.name=="Tuanzi"
    @test c.seniority==1
end
