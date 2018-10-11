```@meta
CurrentModule=Hamiltonian.Utilities.NamedVector
```

# NamedVector

A named vector is similiar to a named tuple, which associate each of its values with a name, represented with a `Symbol`. Unlike named tuples, the values of a named vector can be modified. Yet the names of a named vector cannot be changed. With experiences of `namedtuple` in Python, we treat the names of a named vector as type traits, by overloading the `Base.filednames` function, so that all instances of a certain concrete named vector share the same names, which is different from the predefined `NamedTuple` in Julia. For users who want to define their own concrete named vector, we recommend the use of the macro `@namedvector`.

```@autodocs
Modules=[NamedVector]
Order=  [:module,:constant,:type,:macro,:function]
```
