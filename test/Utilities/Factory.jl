using Test
using Hamiltonian.Utilities.Factory

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

@testset "Parameter" begin
    parameter=Parameter(:Int)
    @test parameter.name==:Int
    @test parameter.type==:nothing
    @test parameter()==:Int

    parameter=Parameter(:(T<:Real))
    @test parameter.name==:T
    @test parameter.type==:Real
    @test parameter()==:(T<:Real)

    parameter=Parameter(:(<:Real))
    @test parameter.name==:nothing
    @test parameter.type==:Real
    @test parameter()==:(<:Real)
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

@testset "Block" begin
    block=@rmlines Block(:(x=1;y=2))
    @test block.body==Any[:(x=1),:(y=2)]
    @test block()==Expr(:block,:(x=1),:(y=2))

    block=@block a=1
    @push! block b=2 c=3 d=4
    @test block.body==[:(a=1),:(b=2),:(c=3),:(d=4)]
end

@testset "FunctionFactory" begin
    ff=FunctionFactory(:(fx(x::Int,y::Int;choice::Function=sum)=choice(x,y)))
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::Int),:(y::Int)])
    @test ff.kwargs==Argument.([:(choice::Function=sum)])
    @test ff.rtype==:Any
    @test ff.params==Parameter[]
    @test ff.body==Block(:(begin choice(x,y) end))

    ff=FunctionFactory(:(
        function fx(x::T,y::T;choice::Function=sum) where T
            choice(x,y)
        end
    ))
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::T),:(y::T)])
    @test ff.kwargs==Argument.([:(choice::Function=sum)])
    @test ff.rtype==:Any
    @test ff.params==Parameter.([:T])
    @test ff.body|>rmlines==Block(:(begin choice(x,y) end))

    ff=@functionfactory (fx()::T)=nothing
    @addargs! ff x::T y::T z::T
    @addkwargs! ff sign::Int=1 choice::Function=sum
    @addparams! ff T<:Number
    @extendbody! ff result=choice(x,y,z) result*=sign
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::T),:(y::T),:(z::T)])
    @test ff.kwargs==Argument.([:(sign::Int=1),:(choice::Function=sum)])
    @test ff.rtype==:T
    @test ff.params==Parameter.([:(T<:Number)])
    @test ff.body|>rmlines==Block(:nothing,:(result=choice(x,y,z)),:(result*=sign))
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
    @test tf.params==Parameter.([:T])
    @test tf.supertype==:(Parent{T} where T <: Real)
    @test tf.fields==Field.([:(address::String),:(phone::Vector{Int}),:(properties::Vector{T})])
    @test tf.constructors==FunctionFactory[]

    tf=@typefactory struct Child <: Parent{T} where T<:Real end
    @addparams! tf T
    @addfields! tf address::String phone::Vector{Int} properties::Vector{T}
    @addconstructors! tf Child(address::String,phone::Vector{Int},properties::Vector{T}) where T=new(address,phone,properties)
    @test tf.name==:Child
    @test tf.mutable==false
    @test tf.params==Parameter.([:T])
    @test tf.supertype==:(Parent{T} where T <: Real)
    @test tf.fields==Field.([:(address::String),:(phone::Vector{Int}),:(properties::Vector{T})])
end

eval((@typefactory struct Child
        name::String
        seniority::Int
        Child(name::String,seniority::Int=1)=new(name,seniority)
    end)())

@testset "call" begin
    c=Child("Tuanzi")
    @test c.name=="Tuanzi"
    @test c.seniority==1
end
