module TypeTraits

using ..Prerequisites: atol,rtol

export efficientoperations
export MemoryOrder,FOrder,COrder
export forder,corder,indtosub,subtoind

struct EfficientOperations end
"""
    efficientoperations

Indicate that the efficient operations, i.e. "=="/"isequal", "<"/"isless" or "replace", will be used.
"""
const efficientoperations=EfficientOperations()

"""
    ==(::EfficientOperations,o1,o2) -> Bool
    isequal(::EfficientOperations,o1,o2) -> Bool

Compare two objects and judge whether they are eqaul to each other.
"""
@generated function Base.:(==)(::EfficientOperations,o1,o2)
    fcount=o1|>fieldcount
    if fcount==o2|>fieldcount
        if fcount==0
            return :(true)
        else
            expr=:(getfield(o1,1)==getfield(o2,1))
            for i=2:fcount
                expr=Expr(:&&,expr,:(getfield(o1,$i)==getfield(o2,$i)))
            end
            return expr
        end
    else
        return :(false)
    end
end
@generated function Base.isequal(::EfficientOperations,o1,o2)
    fcount=o1|>fieldcount
    if fcount==o2|>fieldcount
        if fcount==0
            return :(true)
        else
            expr=:(isequal(getfield(o1,1),getfield(o2,1)))
            for i=2:fcount
                expr=Expr(:&&,expr,:(isequal(getfield(o1,$i),getfield(o2,$i))))
            end
            return expr
        end
    else
        return :(false)
    end
end

"""
    <(::EfficientOperations,o1,o2) -> Bool
    isless(::EfficientOperations,o1,o2) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@generated function Base.:<(::EfficientOperations,o1,o2)
    n1,n2=o1|>fieldcount,o2|>fieldcount
    N=min(n1,n2)
    expr=n1<n2 ? Expr(:if,:(getfield(o1,$N)==getfield(o2,$N)),true,false) : false
    expr=Expr(:if,:(getfield(o1,$N)<getfield(o2,$N)),true,expr)
    for i in range(N-1,stop=1,step=-1)
        expr=Expr(:if,:(getfield(o1,$i)>getfield(o2,$i)),false,expr)
        expr=Expr(:if,:(getfield(o1,$i)<getfield(o2,$i)),true,expr)
    end
    return expr
end
@generated function Base.isless(::EfficientOperations,o1,o2)
    n1,n2=o1|>fieldcount,o2|>fieldcount
    N=min(n1,n2)
    expr=n1<n2 ? Expr(:if,:(getfield(o1,$N)==getfield(o2,$N)),true,false) : false
    expr=Expr(:if,:(isless(getfield(o1,$N),getfield(o2,$N))),true,expr)
    for i in range(N-1,stop=1,step=-1)
        expr=Expr(:if,:(isless(getfield(o2,$i),getfield(o1,$i))),false,expr)
        expr=Expr(:if,:(isless(getfield(o1,$i),getfield(o2,$i))),true,expr)
    end
    return expr
end

"""
    isapprox(::EfficientOperations,::Val{Names},o1,o2;atol::Real=atol,rtol::Real=rtol) where Names -> Bool

Compare two objects and judge whether they are inexactly equivalent to each other.
"""
@generated function Base.isapprox(::EfficientOperations,::Val{Names},o1,o2;atol::Real=atol,rtol::Real=rtol) where Names
    fcount=o1|>fieldcount
    if fcount==o2|>fieldcount
        if fcount==0
            return :(true)
        else
            exprs=[(name=QuoteNode(name);:(isapprox(getfield(o1,$name),getfield(o2,$name);atol=atol,rtol=rtol))) for name in Names]
            for name in fieldnames(o1)
                name ∉ Names && (name=QuoteNode(name);push!(exprs,:(isequal(getfield(o1,$name),getfield(o2,$name)))))
            end
            return Expr(:&&,exprs...)
        end
    else
        return :(false)
    end
end

"""
    replace(::EfficientOperations,o;kwargs...) -> typeof(o)

Return a copy of the input object with some of the field values replaced by the keyword arguments.
"""
@generated function Base.replace(::EfficientOperations,o;kwargs...)
    exprs=[:(get(kwargs,$name,getfield(o,$name))) for name in QuoteNode.(o|>fieldnames)]
    return :(typeof(o).name.wrapper($(exprs...)))
end

abstract type MemoryOrder end
struct FOrder <: MemoryOrder end
struct COrder <: MemoryOrder end
"""
    forder

Indicate that the convertion between Cartesian index and linear index is using the Fortran order.
"""
const forder=FOrder()
"""
    corder

Indicate that the convertion between Cartesian index and linear index is using the C/C++ order.
"""
const corder=COrder()

@generated tail(ts::NTuple{N}) where N=Expr(:tuple,[:(ts[$i]) for i=2:N]...)
@generated head(ts::NTuple{N}) where N=Expr(:tuple,[:(ts[$i]) for i=1:N-1]...)

"""
    indtosub(dims::Tuple,ind::Int,order::FOrder) -> Tuple
    indtosub(dims::Tuple,ind::Int,order::COrder) -> Tuple

Convert an linear index to Cartesian index. Fortran-order or C-order can be assigned.
"""
function indtosub(dims::Tuple,ind::Int,order::FOrder)
    length(dims)==0 && return ()
    ((ind-1)%dims[1]+1,indtosub(tail(dims),(ind-1)÷dims[1]+1,order)...)
end
function indtosub(dims::Tuple,ind::Int,order::COrder)
    length(dims)==0 && return ()
    (indtosub(head(dims),(ind-1)÷dims[end]+1,order)...,(ind-1)%dims[end]+1)
end

"""
    subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N -> Int
    subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N -> Int

Convert an Cartesian index to linear index. Fortran-order or C-order can be assigned.
"""
function subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::FOrder) where N
    length(dims)==0 && return 1
    (subtoind(tail(dims),tail(inds),order)-1)*dims[1]+inds[1]
end
function subtoind(dims::NTuple{N,Int},inds::NTuple{N,Int},order::COrder) where N
    length(dims)==0 && return 1
    (subtoind(head(dims),head(inds),order)-1)*dims[end]+inds[end]
end

end # module
