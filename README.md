# PyDyCoil

## Introduction
PyDycoil is a Python-based analysis package for quantifying the local flexibility and dynamics of filamentous proteins (here, coiled coil structures) from high-speed atomic force microscopy (HS-AFM) images.

Briefly, the following analyses are implemented in the jupyter notebooks:
- *Local flexibility*: This involves aligning individual coiled coils, normalizing their contour lengths, computing persistence length via a 2D worm-like chain (WLC) fit, and averaging the results over time. An example implementation for the test dataset is provided in the local_lp.ipynb Jupyter notebook.

- *Dynamics*: This analysis involves tracking time-resolved fluctuations in coiled coils’ end-to-end distance (R) and bending angle (θ), computing mean-square displacements, and fitting them with an equation for a homogeneous semiflexible filament to extract effective persistence length, internal friction, and relaxation times, followed by averaging and rescaling across conformations. An example implementation for the test dataset is provided in the MSD.ipynb Jupyter notebook.

An example dataset of interpolated X and Y coordinates from the coiled coil trajectory (S configuration) of the bacterial MR complex is provided in the S/ folder.
Each .csv file contains the interpolated X and Y coordinates for each time point (frame) in the trajectory.

## Step-by-step guide
- Clone the repository
```
git clone https://github.com/DyNaMo-INSERM/PyDyCoil.git
cd PyDyCoil
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

- Launch Jupyter Notebook 
```
jupyter notebook
```
- Open and run the notebooks. The repository contains the example dataset and the following analysis notebooks:
 `local_lp.ipynb` – Local flexibility analysis  
 `MSD.ipynb` – Dynamics analysis

### Preprint: 
DOI: 10.64898/2026.07.21.739852
## Citation
If you use PyDyCoil in your research, please cite this repository.
## Contact
If you have any suggestions, comments or encounter any issues, please write to us: yogesh.saravanan@inserm.fr, prithwidip@gmail.com or felix.rico@inserm.fr 



