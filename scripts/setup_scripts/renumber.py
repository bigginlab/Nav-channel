import MDAnalysis as mda

# ref_file = "combined.gro"
ref_file = "0_protein_only.pdb"

# Load in the topology from tleap output files 

u = mda.Universe(ref_file)
dimensions = u.dimensions

# Nav
chain1 = u.select_atoms("index 0-4960")
new_chain1_resids = [i for i in range(118, 430 + 1)]


print(chain1.residues.resids)
print(new_chain1_resids)
print(chain1.residues.resnames)

chain1.residues.resids = new_chain1_resids

chain2 = u.select_atoms("index 4961-8966")
new_chain2_resids = [i for i in range(698, 945 + 1)]
chain2.residues.resids = new_chain2_resids

chain3 = u.select_atoms("index 8967-18697")
new_chain3_resids = [i for i in range(1187, 1782 + 1)]
chain3.residues.resids = new_chain3_resids

# # Lipids
# lipids = u.select_atoms("resname POPC")
# new_lipid_resids = [i for i in range(1, len(lipids.residues.resids) + 1)]
# lipids.residues.resids = new_lipid_resids

# Solvent
#solvent = u.select_atoms("resname Na+ or resname Cl-")
#new_solvent_ion_resids = [i for i in range(1, len(solvent.residues.resids) + 1)]
#solvent.residues.resids = new_solvent_ion_resids

# Create the new system by merging each universe
#new_system = mda.Merge(rbd, rbd_glycans, s309_a, s309_b, solvent)

new_system = mda.Merge(chain1, chain2, chain3)#, lipids)

# Name each chain
new_system.segments.segids = ['A', 'B', 'C']#, 'D']

new_system.dimensions = dimensions

# Write out the new system
#new_system.atoms.write("renumbered.gro")
new_system.atoms.write("renumbered.pdb")
