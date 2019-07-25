```@meta
CurrentModule=QuantumLattices.Interfaces
```

# Interfaces

This module contains the generic functions that are extended by the package.

Due to the multi-dispatch feature of Julia, generic functions can be extended by local methods for different types. However, a local definition of a method also claims a new generic function if the generic function is not imported to the current scope, thus ruins the definitions in other modules. Therefore, it is quite necessary to predefine the common generic functions in a separate module, so that other modules can extend them with their own by a simple import.

## Manual

```@autodocs
Modules = [Interfaces]
Order = [:module, :constant, :type, :macro, :function]
```
