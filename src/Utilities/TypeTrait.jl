module TypeTrait

export efficientoperations

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
    replace(::EfficientOperations,o;kwargs...) -> typeof(o)

Return a copy of the input object with some of the field values replaced by the keyword arguments.
"""
@generated function Base.replace(::EfficientOperations,o;kwargs...)
    exprs=[:(get(kwargs,$name,getfield(o,$name))) for name in QuoteNode.(o|>fieldnames)]
    return :(typeof(o).name.wrapper($(exprs...)))
end

end # module
