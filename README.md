# elmer_mesh_maker

This module makes 2D meshes of different resolutions for ElmerIce, with customised labeled boundaries. ElmerIce will extrude the mesh into 3D.

**See jupyter notebook tutorial.**

Process summary:
1. msh = mesh()
2. msh.load_perimeter(perimeter_filename)
3. msh.set_resolutions(list_of_resolutions)
4. msh.make_perimeters_newres()
5. msh.make_meshes()
6. (optional) msh.export()  
7. msh.print_elmergrid_command()
8. msh.set_corners(list_of_coordinates)
9. msh.label_boundaries(list_of_labels)
