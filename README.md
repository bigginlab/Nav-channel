# nav-alt-splice
Collection of scripts for the setup and analysis of simulations relating to the publication: [Spliced isoforms of the cardiac Nav1.5 channel modify channel activation by distinct structural mechanisms](https://rupress.org/jgp/article-abstract/154/5/e202112906/213074/Spliced-isoforms-of-the-cardiac-Nav1-5-channel?redirectedFrom=fulltext)

# Directory structure
## analysis_scripts
This contains scripts associcated with the geenral analysis of simulations. Typically these produce `.pkl` files that can be used by plotting scripts located in the `plotting_scripts` directory.

An example of `calc_distances.py`:

```
python scripts/analysis_scripts/calc_distances.py -pdb_file input_file.gro -traj_path path/to/trajectories/ -residue_pair 219 208 -plot_hist True -cog True -system_name Double -save_name min_dist_R219_E208
```

## plotting_scripts
This contains a collection of scripts to plot all results together.

An example of `plot_all_distancs.py`:

```
python scripts/plotting_scripts/plot_all_distances.py -data_files 0_WT/207_208_dist_H_carboxy.pkl 1_D211K/207_208_dist_H_carboxy.pkl 2_T207S/207_208_dist_H_carboxy.pkl 3_double_D211K_T207S/207_208_dist_H_carboxy.pkl -save_name min_207_208_H_carboxy_distance_v2 -xlabel "T/S207  - E208 minimum distance" -kde True
```

## pymol_scripts
This contains a collection of PyMol scripts used to generate publication-ready figures of key residues in the voltage-sensor domain.

## setup_scripts
This contains a helper script to remove water present in the membrane when creating an membrane-embedded protein system.