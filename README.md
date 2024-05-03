# Post-ablation image analysis

[![DOI](https://zenodo.org/badge/791883805.svg)](https://zenodo.org/doi/10.5281/zenodo.11068961)

Computes optical flow statistics for ablation data in microscopy images

The position of the ablation available computed from a screenshot (the name of the file starts with AP or DV indicating the orientation of the cut).
The cut is identified and a region of interest is defined around the cut.

<img src=src/Cut.png width="300" height="300"><img src=src/ROI.png width="300" height="300">

Optical flow is computed (Lukas-Kanade algorithm or Brox). The component orhogonal to the cut is calculated.
Mean values inside the ROI are computed (for all and for bright pixels) and  written in a file called Summary.csv in the folder \dara\Result.

<img src=src/OF.png width="300" height="300"><img src=src/quiver2.png width="300" height="300">

This code was used for:

Guy B. Blanchard, Elena Scarpa, Leila Muresan, Bénédicte Sanson
Mechanical stress combines with planar polarised patterning during metaphase to orient embryonic epithelial cell divisions
doi: https://doi.org/10.1101/2023.07.12.548728

The code was tested with Matlab 2022a. An example dataset is provided to be run with the script OFAnalysis.m . 

### external dependency
The source code provided in the folder extrn belongs to Visesh Chari.
Please cite the following publication:
Thomas Brox, Andres Bruhn, Nils Papenberg, Joachim Weickert,
High accuracy optical flow estimation based on a theory for warping,
T. Pajdla and J. Matas (Eds.), European Conference on Computer Vision (ECCV),
Springer, LNCS, Vol. 3024,  25-36, May 2004. 
