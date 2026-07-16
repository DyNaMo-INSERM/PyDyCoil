# PyDyCoil

## Introduction
PyDycoil is a Python-based analysis package for quantifying the local flexibility and dynamics of filamentous proteins (here, coiled-coil structures) from high-speed atomic force microscopy (HS-AFM) imaging data.

Briefly, the following analyses are implemented in the jupyter notebooks:
- *Local flexibility*: This involves aligning individual coiled coils, normalizing their contour lengths, computing persistence length via a 2D worm-like chain (WLC) fit, and averaging the results over time. An example implementation for the test dataset is provided in the local_lp.ipynb Jupyter notebook.

- *Dynamics*: This analysis involves tracking time-resolved fluctuations in coiled coils’ end-to-end distance (R) and bending angle (θ), computing mean-square displacements, and fitting them with an equation for a homogeneous semiflexible filament to extract effective persistence length, internal friction, and relaxation times, followed by averaging and rescaling across conformations. An example implementation for the test dataset is provided in the MSD.ipynb Jupyter notebook.

An example dataset of interpolated X and Y coordinates from the coiled-coil trajectory (S configuration) of the bacterial MR complex is provided in the S/ folder.
Each .csv file contains the interpolated X and Y coordinates for each time point (frame) in the trajectory.

## Step-by-step guide
- Clone the repository
```
git clone https://github.com/DyNaMo-INSERM/PyDyCoil.git
```
- Create an environment with python 3.9
```
conda create -n yourenvname python=3.9 
conda activate yourenvname
```
- Install the dependencies from requirements.txt
```
pip install -r requirements.txt
```

### How to use ⚠️

1. Local flexibility analysis: Use the script *constants.py* together with the main analysis script *XXX.py*. Before running main script, update the file paths in the script to point to the directory where the example CSV files are stored.

2. Dynamics analysis: In the script 'fit_MSD_tau_new', the coiled coil trajectories uses a WLC model. It loads the experimentally measured MSD data from multiple coiled coil configurations and groups trajectories based on individual labelling. For each configuration, the script fits the MSD as a function of time to a theoretical semiflexible polymer relaxation model in which the relaxation modes follow $\tau_n = \tau_1/n^4$, to extrcat the characteristic relaxation time ($\tau$) and persistence length ($l_p$). Using these fitted parameters together with the measured contour length ($L$), it calculates the effective friction coefficient ($\zeta$) and the characteristic fluctuation amplitude given by $r_c = \frac{L^4}{90l_p^2}$. The script then rescales the raw MSD ($t/\tau$ and $\Delta R/r_c$) for direct comparison across datasets (universal curve). ~~In addition, it incorporates independently calculated ($l_p$) and curvature measurements to compute structural descriptors such as mean, median, minimum, maximum, and standard deviation for each filament~~. Finally, the script ~~generates summary statistics, compares the modified $\tau$-based fitting approach with the original $\zeta$-based method, applies the same workflow to human mutant datasets,~~ produces final tables for downstream analysis.


## Acknowledgements
This project has received funding from the Human Frontier Science Program (grant no. RGP0056/2018) and the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program (grant no. 772257).
## Publications ⚠️
Further information about data processing and analysis can be found here:
### Preprint: 
## Citation
If you use PyDyCoil in your research, please cite this repository.
## Contact
If you have any suggestions, comments or encounter any issues, please write to us: yogesh.saravanan@inserm.fr, prithwidip@gmail.com or felix.rico@inserm.fr 



