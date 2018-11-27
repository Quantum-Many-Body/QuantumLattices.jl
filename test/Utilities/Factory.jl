using Hamiltonian.Utilities.Factory

@testset "Inference" begin
    @test (@inference Vector{Tuple{Int,String}})(UnEscaped(:Vector,:Tuple,:Int,:String))==:(Vector{Tuple{Int,String}})
    @test (@inference Type{<:Real})(UnEscaped(:Type,:Real))==:(Type{<:Real})
    @test Inference(:(SVector{N,<:Real}))|>string=="Inference(\n  head:   curly\n  name:   SVector\n  params: Inference[N, <:Real]\n)"
    @test replace(@inference(Vector{<:Real}),name=:AbstractVector)==@inference(AbstractVector{<:Real})
end

@testset "Argument" begin
    argument=Argument(:(phone::Vector{Int}=[9,1,1]))
    @test argument.name==:phone
    @test argument.type==Inference(:(Vector{Int}))
    @test argument.slurp==false
    @test argument.default==:([9,1,1])
    @test argument(MixEscaped(UnEscaped(:Vector,:Int)))==Expr(:kw,:(phone::Vector{Int}),:[9,1,1])

    argument=Argument(:(properties::Int...))
    @test argument.name==:properties
    @test argument.type==Inference(:Int)
    @test argument.slurp==true
    @test argument.default==nothing
    @test argument(MixEscaped(UnEscaped(:Int)))==:(properties::Int...)

    argument1=@argument address::String
    argument2=Argument(:(address::String))
    @test argument1==argument2
end

@testset "Parameter" begin
    parameter=Parameter(:T)
    @test parameter.name==:T
    @test parameter.type==nothing
    @test parameter(UnEscaped(:T))==:T

    parameter=Parameter(:(T<:Real))
    @test parameter.name==:T
    @test parameter.type==Inference(:Real)
    @test parameter(UnEscaped(:Real))==:(T<:Real)

    parameter=Parameter(:(<:Real))
    @test parameter.name==nothing
    @test parameter.type==Inference(:Real)
    @test parameter(UnEscaped(:Real))==:(<:Real)

    parameter=@parameter ::Int
    @test parameter.name==nothing
    @test parameter.type==Inference(:Int)
    @test parameter(UnEscaped(:Int))==:(<:Int)
end

@testset "Field" begin
    field=Field(:xs)
    @test field.name==:xs
    @test field.type==Inference(:Any)
    @test field(UnEscaped(:Any))==:(xs::Any)

    field=Field(:(xs::Vector{Int}))
    @test field.name==:xs
    @test field.type==Inference(:(Vector{Int}))
    @test field(UnEscaped(:Vector,:Int))==:(xs::Vector{Int})

    field1=@field xs::Vector{Int}
    field2=Field(:(xs::Vector{Int}))
    @test field1==field2
end

@testset "Block" begin
    block=@rmlines Block(:(x=1;y=2))
    @test block.body==Any[:(x=1),:(y=2)]
    @test block(Escaped())==Expr(:block,:(x=1),:(y=2))

    block=@block a=1
    @push! block b=2 c=3 d=4
    @test block.body==[:(a=1),:(b=2),:(c=3),:(d=4)]
end

@testset "FunctionFactory" begin
    ff=FunctionFactory(:(fx(x::Int,y::Int;choice::Function=sum)=choice(x,y)))
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::Int),:(y::Int)])
    @test ff.kwargs==Argument.([:(choice::Function=sum)])
    @test ff.rtype==Inference(:Any)
    @test ff.whereparams==Parameter[]
    @test ff.body==Block(:(begin choice(x,y) end))

    ff=FunctionFactory(:(
        function fx(x::T,y::T;choice::Function=sum) where T
            choice(x,y)
        end
    ))
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::T),:(y::T)])
    @test ff.kwargs==Argument.([:(choice::Function=sum)])
    @test ff.rtype==Inference(:Any)
    @test ff.whereparams==Parameter.([:T])
    @test ff.body|>rmlines==Block(:(begin choice(x,y) end))

    ff=@functionfactory (fx()::T)=nothing
    @addargs! ff x::T y::T z::T
    @addkwargs! ff sign::Int=1 choice::Function=sum
    @addwhereparams! ff T<:Number
    @extendbody! ff result=choice(x,y,z) result*=sign
    @test ff.name==:fx
    @test ff.args==Argument.([:(x::T),:(y::T),:(z::T)])
    @test ff.kwargs==Argument.([:(sign::Int=1),:(choice::Function=sum)])
    @test ff.rtype==Inference(:T)
    @test ff.whereparams==Parameter.([:(T<:Number)])
    @test ff.body|>rmlines==Block(:nothing,:(result=choice(x,y,z)),:(result*=sign))
end

@testset "TypeFactory" begin
    tf=TypeFactory(:(
        struct Child{T<:Real} <: Parent{T}
            address::String
            phone::Vector{Int}
            properties::Vector{T}
        end
    ))
    @test tf.name==:Child
    @test tf.mutable==false
    @test tf.params==Parameter.([:(T<:Real)])
    @test tf.supertype==Inference(:(Parent{T}))
    @test tf.fields==Field.([:(address::String),:(phone::Vector{Int}),:(properties::Vector{T})])
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

eval((
    @typefactory struct Child
        name::String
        seniority::Int
        Child(name::String,seniority::Int=1)=new(name,seniority)
    end
)(MixEscaped(UnEscaped(:String,:Int,:Any))))

@testset "call" begin
    c=Child("Tuanzi")
    @test c.name=="Tuanzi"
    @test c.seniority==1
end
