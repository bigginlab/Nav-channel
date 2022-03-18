import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import sys

gro_file = sys.argv[1]

u = mda.Universe(gro_file)

p_ul = []
p_ll = []

sol_to_remove = []

sol_com_list = []

# Find leaflets
L = LeafletFinder(u, "name P")

p_upper = L.groups(0)
p_lower = L.groups(1)

# Update leaflet list with the postions of the P atoms contained within them
for i, p_atom in enumerate(p_upper):

    p_ul.append(p_upper.atoms.positions[i][2])

for i, p_atom in enumerate(p_lower):

    p_ll.append(p_lower.atoms.positions[i][2])

# Get waters
sol_particles = u.select_atoms("name OW")

sol_particles_resids_com_dict = dict.fromkeys(u.select_atoms("resname SOL and name OW"))

# Find the average z coord of each leaflet
ul = sum(p_ul) / float(len(p_ul))
ll = sum(p_ll) / float(len(p_ll))


# Get the z position of the centre of mass of each water molecule
print("Analysing water positions...")
for i, water in enumerate(sol_particles):

    water_z = (sol_particles.positions[i])[2]

    # assign the z coordinate to the O and both H's of the water molecule for extraction later
    sol_particles_resids_com_dict[sol_particles[i]] = water_z

# Extract the resids of the water molecules who's z coordinate of their CoM is not in the bilayer.
print("Removing waters in bilayer, this could take a while...")
sol_resids_not_in_bilayer = [
    k for k, v in sol_particles_resids_com_dict.items() if v <= ll or v >= ul
]

# select everything but water
new_gro_file = u.select_atoms("protein or resname POPC or name NA or name CL")

# make new selection of protein, lipid, ions and water outside of the bilayer
for atom in sol_resids_not_in_bilayer:

    new_gro_file += u.select_atoms("resid " + str(atom.resid))

# write out new selection
new_gro_file.write(str(gro_file.split(".")[0]) + ".gro")

print("New .gro file written!")
