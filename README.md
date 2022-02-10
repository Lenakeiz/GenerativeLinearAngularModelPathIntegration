# Data estimation modelling on the path integration data
Repository modelling of path integration error

Part of my research fellow work in the [Space and Memory lab](https://www.ucl.ac.uk/icn/research/research-groups/space-memory) at UCL

Implemented using **Matlab 2019b** but later versions should work as well

## Respository structure
The main script to execute is DataFitting.m in the root folder
Following folders and subfolder should be added to the Matlab path

- Fitting : contains function used for the data fitting
- Visualization : functions used for debugging the solver and generate modelled points
- Group Comparisons : function to assess differences in fitted parameters from two samples
- Output/Plots : Folders where the graphs generated are saved

## Required packages and dependencies
### Added through the matlab Add-on
- Optimization Toolbox
- Global Optimization Toolbox
- Statistics and Machine Learning Toolbox
- (Parallel Computing Toolbox)
### External packages
- [Circular Statistics Toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
- [notBoxPlot](https://uk.mathworks.com/matlabcentral/fileexchange/26508-notboxplot)

## Short introduction
Path integration is the ability to maintain track of specific position while navigating an environment. The function is known to be affected by systematic errors.

Models of path integration have been developed in the past to outline the origin of these errors. Homing vetor model assume the navigator continuosly updates his postion to keep track
where is it in the environment. By assuming memory decay or leaky integration the errors can be model this way (see [Stangl et al](https://www.nature.com/articles/s41467-020-15805-9)).

Encoding Error models instead assume that the errors are due to errors in encoding the outbound paths. See [Harootonian et al](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007489)
for an example. Both models have been proved to be valid in behavioural studies ([Wiener et al](https://link.springer.com/article/10.1007%2Fs00221-010-2460-7)).
 
Path integration is know to be supported by activity of grid/head direction cells found in the entorhinal cortex, one of the earliest regions affected by AD.

By modelling path integration errors of people at risk of AD we hope to get a better understanding of how the deficits are expressed in navigation and 
ideally that could inform us on an optimal metric to be used for monitoring AD progression or eventually to optimize behavioural tasks.

## Algorithm
The algorithm is an Encoding Error data estimation model.

Data are took from the relative [published article](https://academic.oup.com/brain/article/142/6/1751/5497752) on a path integration task tested on young, older and mildcognitive impairment patients.

Given the problem of solving a triangular path integration problem we model a set of point at the end of each walked segment of the triangle accounting for internal errors in linear and angular estimation.
Our model describes in this way the "mental or thought" space of the participant while he tries to solve the problem.

Errors are modelled using Gaussian distribution for the linear errors and Von Mises distributions for the angular errors that are generated at the end of each segment of the walked path.
The distributions are described in the following image.


![alt text](https://github.com/Lenakeiz/PathIntegrationDataEstimationModelling/blob/main/ReadmeImages/DistributionsAtTheEndOfEachSegment.png "Fitted Distributions at each segment")


The model generates internal space points x<sub>3</sub>' for the real coordinate x<sub>3</sub> described in the following image.


![alt text](https://github.com/Lenakeiz/PathIntegrationDataEstimationModelling/blob/main/ReadmeImages/EstimationForFinalPoint.png "Estimation of final point")


Note that X<sub>3</sub>' represents the participant's guess for the original starting location x<sub>0</sub>, and hence x<sub>3</sub>' is the estimation of x<sub>0</sub>.
Our assumption is that the participant thinks they are back at the original starting location (x<sub>0</sub>) in their internal space, and hence x<sub>3</sub>' = x<sub>0</sub>. However, it is also true that the origin of the paths (x<sub>0</sub>) are the same in the internal and real space (it would not make sense to have an internal error at that point in time).

Therefore to fit the best parameters for our model we maximise the likelihood of the interal datapoints generated at the final point to be close to zero (each path starts at [0,0]).

The likelihood is calculated as the product of datapoints for each trial as we fit parameters per participant.

We used a solver to try to fit the model to our data. There is a section below dedicated to this but at the moment our current preference is the [pareto search](https://uk.mathworks.com/help/gads/paretosearch.html). 

## Running the code
The main script is called DataFitting. You can execute the code blockwise to see the different steps.

# Code structure

## Variables
In the Datafitting.m script there are a set of different global variables that will customize each run. The most important ones are

```
global DEBUG;                 %If 1 fit only one participant
global CATCH_RADIUS;          %Catch radius in meters from the origin that will be used by the fitter to constraint the point close to zero.
global FITTED_PARAMS;         %Should be 3,4 or 5. Number of fitted parameters, useful when switching between different models (no proportional noise, proportional sigma, proportional sigma and beta)
global GENERATE_SEGMENTS;     %Should be 0 or 1. When 1 then generate also points at the end of each walked segment rather than only the final points.
global UNROTATE_MENTAL_SPACE; %Should be 0 or 1. When 1 we will recontruct the real space points from the mental space by basically "undoing" the effect of the fitted parameters on the calculation
```

## Main flow
The code inside DataFitting script has the following macro areas:

1. Loading data from PIT experiment
2. Change the outbund paths to be consistend with having the first segment poiting eastward and then the first turning point to be always anticlockwise. Starting location is always the origin (0,0).
    *This is necessary to make sure all paths are consistent with each other to make calculations easier*. 
3. Fit the paramters
4. Plot results
5. Compare parameters

## Flow Details
The main parameters fitting are performed in **PerformGroupFit** function.

The fuction creates structures for holding the lenght of the segments of each side of the triangle, the angular displacement at the end of each walked length
(first angular displacement is always zero).

For each participant the function then calls the **FitData** function. Within fit data we prepare the solver to minimize the log likelihood estimation. For different solvers see the information below.
Each solver requires different inputs: one parameter is *pointer* to the function that will calculate the log likehood,
other parameters are the initial random estimations of the parameters we are trying to fit in.
The solver will automatically call that function with different parameters estimates. 

The standard estimate function is called *Estimate*. It receives the parameters that needs to be fitted and the trials for each participants. Calculate the estimates for the final positions and calculate the likelihood as a function return.

## Solver
Data model has been fitted using solvers from the Matlab local and global optimization toolbox.

[Fmincon](https://uk.mathworks.com/help/optim/ug/fmincon.html), used widely, is based on gradient descent on the error to look for the best solution.

[Fminsearch](https://uk.mathworks.com/help/matlab/ref/fminsearch.html) is similar but the inner algorithm is different.
The main problem with this solver is the fact are best optimized for convex problems for which we know there is maximum minimum of the function.

In our case there is non linearity between the combination of the parameters and as such a linear solver like fmincon might not be the best option.
The [pareto search](https://uk.mathworks.com/help/gads/paretosearch.html) overcome the non-linear problem by searching into the grid of the parameters space at the cost of being really expensing in terms of computations.

## Plots
PlotGroupData plots the points sampled using the best fitted parameters for each participant. It receives as an input the number of samples we are going to get from the distributions, the participant ID and how many trials for that participant
we are going to plot. Internally it makes also of some of the global variables, so keep that in mind when using it. In particular, 

```CALCULATE_MEAN_TRAJECTORIES``` plot the line corresponding following the mean of the sampled distributions. 

```GENERATE_SEGMENTS``` plot the lines to create the triangles associated with the sampled data points.

If ```UNROTATE_MENTAL_SPACE``` is set to 1, once the data points are sampled the effect of *g* (gain rotation) will be undone. This is used for a sort of sanity check about the model.

## Model Variations
Model with 4 and 5 parameters have been created by allowing the sigma and k to have two contributes on the generation of the distribution: a fixed and a proportional to the walked segment contribute.
For this we allowed sigma first (4 parameters generation) and then we allowed sigma and k (5 parameters generation)
