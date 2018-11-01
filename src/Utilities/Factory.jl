module Factory

export typefactory,functionfactory

struct TypeFactory
    mutable::Bool
    type::Union{Symbol,Expr}
    supertype::Union{Symbol,Expr}
    body::Union{Symbol,Expr}
end


function typefactory

end

function functionfactory

end

end #module
