import meshio
import sys

if len(sys.argv) < 2:
    print("Usage: python msh_AnyTo22.py <mesh_name>")
    sys.exit(1)


mesh_name = sys.argv[1]
if not mesh_name.endswith(".msh"):
    print("Mesh filename must end with .msh extension.")
    print(f"Provided name: {mesh_name}")
    sys.exit(1)

out_mesh_name = f"{mesh_name[:-4]}_gmsh22.msh"

mesh = meshio.read(mesh_name, file_format="gmsh")

# if isinstance(mesh, meshio._mesh.Mesh):
#     print(mesh.cell_data)
#     print(mesh.cell_sets)
#     print(mesh.cells)

#     # mesh.cell_sets_to_data()
#     # print(mesh.cell_data)
#     # mesh.cell_data["gmsh:physical"] = mesh.cell_data["PT_GM-PT_CSF-PT_WM-PT_VENT"]
#     # del mesh.cell_data["PT_GM-PT_CSF-PT_WM-PT_VENT"]
#     # print(mesh.cell_data)

print("Starting write.")
meshio.write(out_mesh_name, mesh, file_format="gmsh22", binary=False)

print(f"Done. Mesh written to {out_mesh_name}")