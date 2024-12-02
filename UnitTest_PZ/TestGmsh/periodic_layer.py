import gmsh
import math
import os
import sys

gmsh.initialize()

gmsh.model.add("periodic_test")
x = y = z = 0
dx = 0.5
dy = 1
rect_tag = gmsh.model.occ.add_rectangle(x,y,z,dx,dy)
gmsh.model.occ.synchronize()

elsize = 0.5*dy
gmsh.model.mesh.set_size(gmsh.model.get_entities(0), elsize)

dim = 2
bnds = gmsh.model.get_boundary([(dim,rect_tag)],oriented=False)
mclist = [gmsh.model.occ.get_center_of_mass(d,t) for d,t in bnds]
min_x = bnds[[i for i in range(len(mclist)) if mclist[i][0] == 0][0]][1]
max_x = bnds[[i for i in range(len(mclist)) if mclist[i][0] == 0.5][0]][1]

other_bnds = [t for _,t in bnds]
other_bnds.remove(min_x)
other_bnds.remove(max_x)
dim = 1

affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
pos = {"dx": 3, "dy": 7, "dz": 11}


affine[pos["dx"]] = dx
affine[pos["dy"]] = 0
gmsh.model.mesh.set_periodic(dim,[max_x], [min_x],affine)



gmsh.model.add_physical_group(2, [rect_tag], 1, name="surf")
gmsh.model.add_physical_group(1, [min_x], 11, name="min_x")
gmsh.model.add_physical_group(1, [max_x], 12, name="max_x")
gmsh.model.add_physical_group(1, other_bnds, 13, name="other")


gmsh.model.mesh.generate(2)
gmsh.write("periodic_layer.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
