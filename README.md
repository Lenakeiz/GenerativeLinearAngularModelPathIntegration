Generative Linear-Angular Model of Path Integration (GLAMPI)  
======

[Figure 1](https://raw.githubusercontent.com/Lenakeiz/VectorAdditionModel/blob/main/Images/Fig1a.png)

## Description 
This repository contains the full dataset and the scripts for generating the analysis and the figures from _"Dissociating linear and angular path integration in ageing and Alzheimer’s disease"_[^1]  
The anonymised dataset for healthy elderly and MCI patients is the one from the _"Differentiation of mild cognitive impairment using an entorhinal cortex-based test of virtual reality navigation"_[^2]

## Installation
Clone the repository anywhere in your machine using `git clone`. 
Open the folder in Matlab by making sure to add folders and subfolders to the path.

## Package dependency
Developed using Matlab R2022a.

Requires the following MATLAB toolboxes:

- [Statistics and machine learning](https://uk.mathworks.com/products/statistics.html)
- [Optimization toolbox](https://uk.mathworks.com/products/optimization.html)
- [Global Optimization toolbox](https://uk.mathworks.com/products/global-optimization.html?s_tid=srchtitle_Global%20Optimization%20Toolbox_1)
- [Ecoometric toolbox](https://uk.mathworks.com/products/econometrics.html?s_tid=srchtitle_econometrics%20toolbox_1)
- [Mapping toolbox](https://uk.mathworks.com/products/mapping.html)

We also used a customised version of [Rob Campbell (2023) sigstar](https://github.com/raacampbell/sigstar) which is included already in the project files.

## Use
Scripts to reproduce the plots and analysis in the paper can be found into the **Main** folder. 
Each of the scripts is standalone, meaning it will produce one or more plots from any given figure of the main text/supplementary. 
Plots will not be displayed within Matlab but output will be saved inside **Output** folder.
Some scripts will take longer than others to run (like the one that runs all different types of model). If you want to run the entire analysis, you can do it with `RunAllAnalysis`. 

Below you can find some more information about the analysis process with some references to relevant scripts/functions.

### Configuration

To allow the model to be run with different configurations we have used a `config` struct variable. To see information about properties set in the struct please refer to `GLAMPI_PrepareBaseConfig`. To run a specific version of our generative model the variables `ModelName`,`ParamName` and `NumParams` have to be correctly set. Please see `Fig_2_S3_S4_ModelSelection.m` script to see all available combinations of the model.

### Preprocessing 

Details for the preprocessing of the data is described in the online methods. The script responsible is `GLAMPI_PreprocessData`. Other important functions are `TransformPaths` aligning the paths, `CalculateTrackingPath` for using tracking data to extract information about the participant behaviour while performing the outbound path and `CalculateOoB` to recover relevant information from the trials marked as out of boundary.

### Fitting the model/Behavioural analysis

The entry point for calculating behavioural performances and fitting the model is `GLAMPI`. In `PerformGroupFit` behavioural data is calculated. `FitData` executes the model fitting using the specified model (see `config` struct variable) with the limits for each of the fitted parameters. The fitter use `fmincon` Matlab function, but adopts a global search approach to avoid the local minima problem.

### Visualization
For visualizing trials you can use the `VisualizeRealtimeTrackingData` function.

---
[^1]: Andrea Castegnaro *, Zilong Ji *, Katarzyna Rudzka, Dennis Chan, Neil Burgess _( * ) Equal contributions_

[^2]: David Howett *, Andrea Castegnaro *, Katarzyna Krzywicka, Johanna Hagman, Deepti Marchment, Richard Henson, Miguel Rio, John A King, Neil Burgess, Dennis Chan _( * ) Equal contributions_