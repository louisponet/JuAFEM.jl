# # Cantilever Beam
# ```@setup @__NAME__
# cp(joinpath("@__DIR__", "..", "..", "..", "examples", "figures", "cant.png"), "cantilever.png")
# ```
# ![](cantilever.png)

# ## Introduction

# The heat equation is the "Hello, world!" equation of finite elements.
# Here we solve the equation on a unit square, with a uniform internal source.
# The strong form of the (linear) heat equation is given by
# ```math
# - \nabla \cdot (k \nabla u) = f  \quad x \in \Omega,
# ```
# where $u$ is the unknown temperature field, $k$ the heat conductivity,
# $f$ the heat source and $\Omega$ the domain. For simplicity we set $f = 1$
# and $k = 1$. We will consider homogeneous Dirichlet boundary conditions such that
# ```math
# u(x) = 0 \quad x \in \partial \Omega,
# ```
# where $\partial \Omega$ denotes the boundary of $\Omega$.

# The resulting weak form is given by
# ```math
# \int_{\Omega} \nabla v \cdot \nabla u \ d\Omega = \int_{\Omega} v \ d\Omega,
# ```
# where $v$ is a suitable test function.

# ## Commented program

# Now we solve the problem in JuAFEM. What follows is a program spliced with comments.
# The full program, without comments, can be found in the next [section](@ref @__NAME__-plain-program).

# First we will load some packages.
using JuAFEM
using Tensors
using TimerOutputs
using UnicodePlots
const to = TimerOutput()

# ### Grid
# Next we make grid
corner1 = Vec{dim}((0.0, 0.0, 0.0))
corner2 = Vec{dim}((10.0, 1.0, 1.0))
grid = generate_grid(geoshape, (60, 6, 6), corner1, corner2)

# ### FEValues
interpolation_space = Lagrange{dim, refshape, 1}()
quadrature_rule = QuadratureRule{dim, refshape}(order)
cellvalues = CellVectorValues(quadrature_rule, interpolation_space);
facevalues = FaceVectorValues(QuadratureRule{dim-1, refshape}(order), interpolation_space)

# ### DofHandler
dh = DofHandler(grid)
push!(dh, :u, dim) # Add a displacement field
close!(dh)

# Sparsity, assemble only upper half since it is symmetric
@time K = create_symmetric_sparsity_pattern(dh)
fill!(K.data.nzval, 1.0);
spy(K.data)

# ### Boundary conditions
dbcs = ConstraintHandler(dh)
# Add a homogenoush boundary condition on the "clamped" edge
dbc = Dirichlet(:u, getfaceset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], collect(1:dim))
add!(dbcs, dbc)
close!(dbcs)
t = 0.0
update!(dbcs, t)

# ### Material modeling
# Create the stiffness tensor
E = 200e9
ν = 0.3
λ = E*ν / ((1+ν) * (1 - 2ν))
μ = E / (2(1+ν))
δ(i,j) = i == j ? 1.0 : 0.0
g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
C = SymmetricTensor{4, dim}(g)

# ### Assembling
function doassemble{dim}(cellvalues::CellVectorValues{dim}, facevalues::FaceVectorValues{dim},
                         K::Symmetric, grid::Grid, dh::DofHandler, C::SymmetricTensor{4, dim})


    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)

    n_basefuncs = getnbasefunctions(cellvalues)

    fe = zeros(n_basefuncs) # Local force vector
    Ke = Symmetric(zeros(n_basefuncs, n_basefuncs), :U) # Local stiffness mastrix

    t = Vec{3}((0.0, 1e8, 0.0)) # Traction vector
    b = Vec{3}((0.0, 0.0, 0.0)) # Body force
    ɛ = [zero(SymmetricTensor{2, dim}) for i in 1:n_basefuncs]
    @inbounds for (cellcount, cell) in enumerate(CellIterator(dh))
        @timeit to "assem" begin
        fill!(Ke.data, 0)
        fill!(fe, 0)

        reinit!(cellvalues, cell)
        for q_point in 1:getnquadpoints(cellvalues)
            for i in 1:n_basefuncs
                ɛ[i] = symmetric(shape_gradient(cellvalues, q_point, i))
            end
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                δu = shape_value(cellvalues, q_point, i)
                fe[i] += (δu ⋅ b) * dΩ
                ɛC = ɛ[i] ⊡ C
                for j in i:n_basefuncs # assemble only upper half
                    Ke.data[i, j] += (ɛC ⊡ ɛ[j]) * dΩ # can only assign to parent of the Symmetric wrapper
                end
            end
        end

        for face in 1:nfaces(cell)
            if onboundary(cell, face) && (cellcount, face) ∈ getfaceset(grid, "right")
                reinit!(facevalues, cell, face)
                for q_point in 1:getnquadpoints(facevalues)
                    dΓ = getdetJdV(facevalues, q_point)
                    for i in 1:n_basefuncs
                        δu = shape_value(facevalues, q_point, i)
                        fe[i] += (δu ⋅ t) * dΓ
                    end
                end
            end
        end
        global_dofs = celldofs(cell)
        assemble!(assembler, global_dofs, fe, Ke)
        end # timer
    end
    return K, f
end

# ### Solution
reset_timer!(to)
K, f = doassemble(cellvalues, facevalues, K, grid, dh, C);
print_timer(to; linechars = :ascii)

# Modify K and f such that K \ f gives correct boundary conditions
@time apply!(K, f, dbcs)

@time u = cholfact(K) \ f;

# ### Export
# Save file
vtk_grid("cantilever", dh) do vtk
    vtk_point_data(vtkfile, dh, u)
end

# ## [Plain program](@id @__NAME__-plain-program)

# Below follows a version of the program without any comments.
# The file is also available here: [\[@__NAME__.jl\]](@__NAME__.jl)

# @__CODE__
