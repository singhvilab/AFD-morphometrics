# AFD-morphometrics
Quantification of neuron-receptive endings shape and properties from high-resolution imaging data

This repository contains codes used to quantify image datasets of neuron-receptive endings related to the manuscript:
"Glial Parkin inhibits developmental activity-dependent neuron pruning via Rac1 and Jun Kinase", Teets *et al.*, [Singhvi Lab](https://research.fredhutch.org/singhvi/en.html)

## Software requirements

* MATLAB R2023b or higher
* Optimization Toolbox
* Signal Processing Toolbox
* Image Processing Toolbox
* Statistics and Machine Learning Toolbox
* Curve Fitting Toolbox
* Parallel Computing Toolbox

## Dependencies

* [bfmatlab package](https://www.openmicroscopy.org/bio-formats/downloads/) to use Bio-Formats for MATLAB
* [Adaptive Resolution Orientation Space](https://github.com/mkitti/AdaptiveResolutionOrientationSpace) from Mark Kittisopikul

The analysis is launched using the `eteetsWrap()` function.

`minimalBridge` is originally found in the Adaptive Resolution Space repo. This is an updated version with a bug fixed due to function depreciation since R2022.
