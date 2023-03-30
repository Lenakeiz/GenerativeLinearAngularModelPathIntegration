Generative Linear-Angular Model of Path Integration (GLAMPI)  
======

[Figure 1]: https://github.com/Lenakeiz/VectorAdditionModel/blob/main/Images/Fig1a.png "Triangle Completion Task"

## Description 
This repository contains the full dataset and the scripts for generating the analysis and the figures from _"Dissociating linear and angular path integration in ageing and Alzheimer’s disease"_[^1]  
The anonymised dataset for healthy elderly and MCI patients is the one from the _"Differentiation of mild cognitive impairment using an entorhinal cortex-based test of virtual reality navigation"_[^2]

---
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

---
## Use
Scripts to reproduce the plots and analysis in the paper can be found into the **Main** folder. 
Each of the scripts is standalone, meaning it will produce one or more plots from any given figure of the main text/supplementary. 
Plots will not be displayed within Matlab but output will be saved inside **Output** folder.
Some scripts will take longer than others to run (like the one that runs all different types of model).

### Preprocessing 
The preprocessing of the data as described in the 

---
[^1]: Andrea Castegnaro *, Zilong Ji *, Katarzyna Rudzka, Dennis Chan, Neil Burgess _( * ) Equal contributions_

[^2]: David Howett *, Andrea Castegnaro *, Katarzyna Krzywicka, Johanna Hagman, Deepti Marchment, Richard Henson, Miguel Rio, John A King, Neil Burgess, Dennis Chan _( * ) Equal contributions_