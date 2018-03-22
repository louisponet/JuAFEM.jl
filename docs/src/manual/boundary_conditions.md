```@meta
DocTestSetup = :(using JuAFEM)
```

# Boundary Conditions

To every PDE there is boundary conditions. There are different types of boundary
conditions, and they need to be handled in different ways. Below we discuss how to handle
the most common ones, Dirichlet and Neumann boundary conditions, and how to do it `JuAFEM`.

## Dirichlet Boundary Conditions

At a Dirichlet boundary the solution is prescribed to a given value. For the discrete
FE-solution this means that there are some degrees of freedom that are fixed. To be able
to tell which degrees of freedom we should constrain we need the `DofHandler`.

```julia
ch = ConstraintHandler(dh)
```

Check out the following commented examples, which deals with Dirichlet boundary conditions:
 - [Heat Equation](@ref)
 - TODO


## Neumann Boundary Conditions

At the Neumann part of the boundary we know something about the gradient of the solution.

As an example, we could implement this as

```julia
for face in cell
    if onboundary(face, cell) && (cell, face) in getfaceset(grid, "Neumann Boundary")
        for qp in 1:...
            for i in 1:getnbasefunctions()
                fe[i]
            end
        end
    end
end
```

Check out the following commented examples, which deals with Neumann boundary conditions:
 - TODO
