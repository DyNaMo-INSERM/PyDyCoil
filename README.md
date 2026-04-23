# PyDyCoil

## Introduction
PyDycoil is a Python-based analysis package for quantifying the local flexibility and dynamics of filamentous proteins (here, coiled-coil structures) from high-speed atomic force microscopy (HS-AFM) imaging data.
Briefly, the following analyses are implemented in the scripts:
- *Local flexibility*: This involves aligning individual coiled coils, normalizing their contour lengths, computing persistence length via a 2D worm-like chain (WLC) fit, and averaging the results.
- *Dynamics*: This analysis involves tracking time-resolved fluctuations in coiled coils’ end-to-end distance (R) and bending angle (θ), computing mean-square displacements, and fitting them with an equation for a homogeneous semiflexible filament to extract effective persistence length, internal friction, and relaxation times, followed by averaging and rescaling across conformations.

## Step-by-step guide ⚠️
### Installation of python packages⚠️
✅ Ideas

This code was developed using: Python 3.10+  
Required Python libraries:
pandas (reading CSV files), numpy (numerical operations), scipy (curve fitting), matplotlib (plotting), seaborn (statistical visualization)

Recommended installation:

Install Anaconda (that includes both Python and Spyder IDE used in our case) and install the above packages if they are missing.

### How to use ⚠️
All the scripts are avilable in the folder XXX. For local felxibility analysis, open the following scripts XXX

The csv files are provided in the folder XXXX. Currently, chnage the file path to XXX where the csv file exsits on your computer. 

Input CSV files>>>>

need to find out all the csv files as an input.

in the fit_MSD_tau_new, the This script performs a complete analysis pipeline to quantify the dynamic and mechanical properties of filament trajectories using a worm-like chain (WLC) relaxation model. It begins by loading experimentally measured mean squared displacement (MSD) data from multiple filament configurations, including S, C, ring, braid, MRN, monomer, volume variants, and human mutant datasets, and groups trajectories based on individual file and filament labels. For each filament, the script fits the MSD as a function of time to a theoretical semiflexible polymer relaxation model in which the relaxation modes follow τ 
n
​	
 =τ 
1
​	
 /n 
4
 , enabling extraction of the characteristic relaxation time (τ) and persistence length (l 
p
​	
 ). Using these fitted parameters along with the measured contour length (L), it calculates the effective friction coefficient (ζ) and characteristic fluctuation amplitude r 
c
​	
 =L 
4
 /(90l 
p
2
​	
 ). The script then rescales the raw MSD trajectories into dimensionless variables (t/τ and ΔR/r 
c
​	
 ) to enable direct comparison across datasets on a universal basis. In addition, it incorporates independently calculated local persistence length and curvature measurements to compute structural descriptors such as mean, median, minimum, maximum, and standard deviation for each filament. Finally, the pipeline generates summary statistics across experimental conditions, compares the modified τ-based fitting approach with the original ζ-based method, applies the same workflow to human mutant datasets, and produces final tables and plots for downstream analysis, comparison, and publication.




So we can use for example, the interpolated cooridnates of two sets of data from S configuration. 















## Acknowledgements
This project has received funding from the Human Frontier Science Program (grant no. RGP0056/2018) and the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program (grant no. 772257).
## Publications ⚠️
### Preprint: 
## Citation
If you use PyDyCoil in your research, please cite this repository.
## Contact
If you have any suggestions, comments or encounter any issues, please write to us: yogesh.saravanan@inserm.fr, prithwidip@gmail.com or felix.rico@inserm.fr 



