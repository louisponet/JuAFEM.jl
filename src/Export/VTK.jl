cell_to_vtkcell(::Type{Line}) = VTKCellTypes.VTK_LINE
cell_to_vtkcell(::Type{QuadraticLine}) = VTKCellTypes.VTK_QUADRATIC_EDGE

cell_to_vtkcell(::Type{Quadrilateral}) = VTKCellTypes.VTK_QUAD
cell_to_vtkcell(::Type{QuadraticQuadrilateral}) = VTKCellTypes.VTK_BIQUADRATIC_QUAD
cell_to_vtkcell(::Type{Triangle}) = VTKCellTypes.VTK_TRIANGLE
cell_to_vtkcell(::Type{QuadraticTriangle}) = VTKCellTypes.VTK_QUADRATIC_TRIANGLE
cell_to_vtkcell(::Type{Cell{2,8,4}}) = VTKCellTypes.VTK_QUADRATIC_QUAD

cell_to_vtkcell(::Type{Hexahedron}) = VTKCellTypes.VTK_HEXAHEDRON
cell_to_vtkcell(::Type{Tetrahedron}) = VTKCellTypes.VTK_TETRA
cell_to_vtkcell(::Type{QuadraticTetrahedron}) = VTKCellTypes.VTK_QUADRATIC_TETRA

"""
    vtk_grid(filename::AbstractString, grid::Grid; linear=false)

Create a unstructured VTK grid from a `Grid`. Return a `DatasetFile`
which data can be appended to, see `vtk_point_data` and `vtk_cell_data`.
The `linear` keyword can be set to `true` to export only the ``linear''
part of the grid.
"""
function WriteVTK.vtk_grid(filename::AbstractString, grid::Grid{dim,N,T}; linear::Bool=false) where {dim,N,T}
    celltype = cell_to_vtkcell(getcelltype(grid))
    cls = MeshCell[]
    for cell in CellIterator(grid)
        push!(cls, MeshCell(celltype, getnodes(cell)))
    end
    coords = reinterpret(T, getnodes(grid), (dim, getnnodes(grid)))
    return vtk_grid(filename, coords, cls)
end

# function WriteVTK.vtk_point_data(vtk::WriteVTK.DatasetFile, data::Vector{Vec{dim,T}}, name::AbstractString) where {dim,T}
#     npoints = length(data)
#     data = reinterpret(T, data, (dim, npoints))
#     return vtk_point_data(vtk, data, name)
# end

# function vtk_nodeset(vtk::WriteVTK.DatasetFile, grid::Grid{dim}, nodeset::String) where {dim}
#     z = zeros(getnnodes(grid))
#     z[collect(getnodeset(grid, nodeset))] = 1.0
#     vtk_point_data(vtk, z, nodeset)
# end

# function vtk_cellset(vtk::WriteVTK.DatasetFile, grid::Grid{dim}, cellset::String) where {dim}
#     z = zeros(getncells(grid))
#     z[collect(getcellset(grid, cellset))] = 1.0
#     vtk_cell_data(vtk, z, cellset)
# end

module VTK
import WriteVTK

save(args...) = WriteVTK.vtk_save(args...)
grid(args...) = WriteVTK.vtk_grid(args...)
point_data(args...) = WriteVTK.vtk_point_data(args...)
cell_data(args...) = WriteVTK.vtk_cell_data(args...)


end # module VTK



cv.shape_value(qp, i)


shape_value(ci, qp, i)    -> cv.shape_value(qp, i)
shape_gradient(ci, qp, i) -> cv.shape_gradient(qp, i)
function_value(ci, qp, ue)    -> cv.function_value(qp, ue)
getdetJdV(cv, qp)         -> cv.detJdV(qp)

getnquadpoints(cv) -> cv.nquadpoints()
