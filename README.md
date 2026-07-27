# PyDyCoil

## Introduction
PyDycoil is a Python-based analysis package for quantifying the local flexibility and dynamics of filamentous proteins (here, coiled coil structures) from high-speed atomic force microscopy (HS-AFM) videos.

Briefly, the following analyses are implemented in the jupyter notebooks:
- *Local flexibility*: This involves aligning individual coiled coils, normalizing their contour lengths, computing persistence length via a 2D worm-like chain (WLC) fit, and averaging the local persistence length over time. An example implementation for the test dataset is provided in the [local_lp.ipynb](local_lp.ipynb) Jupyter notebook.

- *Dynamics*: This analysis quantifies the dynamics from the fluctuations in the coiled coils end-to-end distance (R) and bending angle (θ), by computing mean-square displacements, and fitting them a theoretical model for a semiflexible filament to extract effective persistence length, internal friction, and relaxation times. The extracted quantities are then averaged, and the dynamics are rescaled across video and different conformations. An example implementation for an example dataset is provided in the [MSD.ipynb](MSD.ipynb) Jupyter notebook.

The example input dataset for these analyses consists of interpolated X and Y coordinates extracted from a video of the coiled coil trajectory (S configuration) of the bacterial MR complex and is provided in the [`S/`](S/) folder. Each `.csv` file in this folder contains the interpolated X and Y coordinates for a single time point (frame) in the trajectory.

## Step-by-Step Guide
- Clone the repository
```
git clone https://github.com/DyNaMo-INSERM/PyDyCoil.git
cd PyDyCoil
```
- Create an environment with python 3.10
```
conda create -n yourenvname python=3.10 
conda activate yourenvname
```
- Install the dependencies from requirements.txt
```
pip install -r requirements.txt
```

- Launch Jupyter Notebook 
```
jupyter lab
```
- Open and run the notebooks. The repository contains the example dataset and the following analysis notebooks:
 [local_lp.ipynb](local_lp.ipynb) – Local flexibility analysis and [MSD.ipynb](MSD.ipynb) – Dynamics analysis

## Citation
If you use **PyDyCoil** in your research, please cite the associated preprint:

Saha P, Saravanan Y, Marchesi A, *et al.* **Scaling the Dynamics of Coiled Coils**. *bioRxiv* (2026).
(https://doi.org/10.64898/2026.07.21.739852)

## Contact
If you have any suggestions, comments or encounter any issues, please write to us: yogesh.saravanan@inserm.fr, prithwidip@gmail.com or felix.rico@inserm.fr 



