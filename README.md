# PyDyCoil

## Introduction
PyDycoil is a Python-based analysis package for quantifying the local flexibility and dynamics of filamentous proteins (here, coiled-coil structures) from high-speed atomic force microscopy (HS-AFM) imaging data.
Briefly, the following analyses are implemented in the scripts:
- *Local flexibility*: This involves aligning individual coiled coils, normalizing their contour lengths, computing persistence length via a 2D worm-like chain (WLC) fit, and averaging the results.
- *Dynamics*: This analysis involves tracking time-resolved fluctuations in coiled coils’ end-to-end distance (R) and bending angle (θ), computing mean-square displacements, and fitting them with an equation for a homogeneous semiflexible filament to extract effective persistence length, internal friction, and relaxation times, followed by averaging and rescaling across conformations.

## Step-by-step guide
### Installation
- The scripts are compatible with Python 3.10 or higher 
- Required Python libraries: pandas (reading CSV files), numpy (numerical operations), scipy (curve fitting), matplotlib (plotting), seaborn (statistical visualization)
- Install Anaconda distribution that includes Python, Spyder IDE (used in our paper) and most libraries preinstalled. If any dependencies are missing in your environment, they can be installed using, for example:
```bash
pip install pandas numpy scipy matplotlib seaborn
```
Download all scripts in the folder named *scripts*, along with the example CSV files provided in the *example* folder, to run the analysis workflow locally on your computer. The example datasets are interpolated coordinates from coiled coil's trajectory from bacterial MR complex.

### How to use ⚠️

1. Local flexibility analysis: Use the script *constants.py* together with the main analysis script *XXX.py*. Before running main script, update the file paths in the script to point to the directory where the example CSV files are stored.

2. Dynamics analysis: In the script 'fit_MSD_tau_new', the coiled coil trajectories uses a WLC model. It loads the experimentally measured MSD data from multiple coiled coil configurations and groups trajectories based on individual labelling. For each configuration, the script fits the MSD as a function of time to a theoretical semiflexible polymer relaxation model in which the relaxation modes follow $\tau_n = \tau_1/n^4$, to extrcat the characteristic relaxation time ($\tau$) and persistence length ($l_p$). Using these fitted parameters together with the measured contour length ($L$), it calculates the effective friction coefficient ($\zeta$) and the characteristic fluctuation amplitude given by $r_c = \frac{L^4}{90l_p^2}$. The script then rescales the raw MSD ($t/\tau$ and $\Delta R/r_c$) for direct comparison across datasets (universal curve). ~~In addition, it incorporates independently calculated ($l_p$) and curvature measurements to compute structural descriptors such as mean, median, minimum, maximum, and standard deviation for each filament. Finally, the script generates summary statistics, compares the modified $\tau$-based fitting approach with the original $\zeta$-based method, applies the same workflow to human mutant datasets, and produces final tables and plots for downstream analysis.


## Acknowledgements
This project has received funding from the Human Frontier Science Program (grant no. RGP0056/2018) and the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program (grant no. 772257).
## Publications ⚠️
### Preprint: 
## Citation
If you use PyDyCoil in your research, please cite this repository.
## Contact
If you have any suggestions, comments or encounter any issues, please write to us: yogesh.saravanan@inserm.fr, prithwidip@gmail.com or felix.rico@inserm.fr 



