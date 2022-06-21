# The magnitude of the seasonal forcing of RSV and implications for vaccination strategies
Fabienne Krauer<sup>1</sup>, Tor Erlend Fjelde<sup>2</sup>, Mihaly Koltai<sup>1</sup>, David Hodgson<sup>1</sup>, Marina Treskova-Schwarzbach<sup>3</sup>, Christine Harvey<sup>4</sup>, Mark Jit<sup>1</sup>, Ole Wichmann<sup>3</sup>, Thomas Harder<sup>3</sup>, Stefan Flasche<sup>1</sup>


<sup>1</sup> [Centre for Mathematical Modelling of Infectious Diseases](https://www.lshtm.ac.uk/research/centres/centre-mathematical-modelling-infectious-diseases), London School of Hygiene & Tropical Medicine, London, UK<br/>
<sup>2</sup> [Computational and Biological Learning, University of Cambridge](http://learning.eng.cam.ac.uk/Public/), UK<br/>
<sup>3</sup> [Robert Koch Institut](https://www.rki.de), Berlin, Germany<br/>
<sup>4</sup> [Health Protection NSW](https://www.health.nsw.gov.au/), NSW Ministry of Health<br/>

## Abstract
**Background**: 
Respiratory syncytial virus (RSV) is a leading cause of respiratory tract infections and bronchiolitis in young children. The seasonal pattern of RSV is shaped by short-lived immunity, seasonally varying contact rates and pathogen viability. The magnitude of each of these parameters is not fully clear. The disruption of the regular seasonality of RSV during the COVID pandemic in 2020 due to control measures, and the ensuing delayed surge in RSV cases provides an opportunity to disentangle these factors and to understand the implication for vaccination strategies. A better understanding of the drivers of RSV seasonality is key for developing future vaccination strategies.<br/> 
**Methods**: 
We developed a mathematical model of RSV transmission, which simulates the sequential re-infection (SEIRRS4) and uses a flexible Von Mises function to model the seasonal forcing. Using MCMC we fit the model to laboratory confirmed RSV data from 2010-2022 from NSW while accounting for the reduced contact rates during the pandemic with Google mobility data. We estimated the baseline transmission rate, its amplitude and shape during RSV season as well as the duration of immunity. The resulting parameter estimates were compared to a fit to pre-pandemic data only, and to a fit with a cosine forcing function. We then simulated the expected shifts in peak timing and amplitude under two vaccination strategies: continuous and seasonal vaccination.<br/> 
**Results**: 
We estimate that RSV dynamics in NSW can be best explained by a high effective baseline transmission rate (2.94/d, 95% CrI 2.72-3.19) and a narrow peak with a maximum 13% increase compared to the baseline transmission rate. We also estimate the duration of post infection temporary but sterilizing immunity to be 412 days (95% CrI 391-434). A cosine forcing resulted in a similar fit and posterior estimates. Excluding the data from the pandemic period in the fit increased parameter correlation and yielded less informative posterior distributions. The continuous vaccination strategy led to more extreme seasonal incidence with a delay in the peak timing and a higher amplitude whereas seasonal vaccination flattened the incidence curves.<br/> 
**Conclusion**: 
Quantifying the parameters that govern RSV seasonality is key in determining potential indirect effects from immunization strategies as those are being rolled out in the next few years.


## Software
The main model code is written in Julia and fitted with [Turing](https://turing.ml/dev/). You can download Julia [here](https://julialang.org/downloads/). 

## Repository content and instructions

This repository contains all Julia and R code for this analysis as well as the input data and the final posterior samples (traces). The main scripts are:

| script name | description |
| :--- | :--- |
| **prep.R** | preparation of input data |
| **results.R** | produces all final figures |
| **scripts/fit.jl** | fits the main model to the weekly RSV cases |
| **scripts/visualize.jl** | to quickly visualize the results of the selected fit(s). Produces a PDF file |
| **scripts/postprocessing.jl** | simulates the posterior predictive of the selected fit(s). Produces several csv files and figures |
| **scripts/vaccsim.jl** | simulates the theoretical disease dynamics under different vaccination scenarios, based on the posteriors from the selected fit(s). Produces several csv files |
| **scripts/loo.jl** | calculated the Pareto-smoothed importance sampling leave-one-out cross-validation (LOO), based on the posteriors from the selected fit(s). |
| **src/models.jl** | Contains the turingmodel |
| **src/differential_equations/** | Contains the ODE models for the fitting (SEIRRS4) and the vaccination simulatino (SEIRRS4Vacc and SEIRRS4Vacc_pulse) |


To get started with any julia file, run the following from the commandline:
```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
```

To start the **fitting process** with the default settings, run the following from the julia console:
```julia
julia> using DrWatson
julia> _args = ["seed", "--betafunc", "--symb", "--verbose"]
julia> include(scriptsdir("fit.jl"))
```

Alternatively you can run the following line in the console:
```
julia --project scripts/fit.jl seed --betafunc=mises --symb --verbose
```

*seed* is any seed no. you prefer, *betafunc* refers to the seasonal forcing function (mises or cosine) and *symb* makes a symbolic version of the init ODE model for faster fitting. *verbose* returns more detailed information of what the script is doing at any time, we recommend to turn it on. `fit.jl` takes multiple arguments to adjust the solver and sampling settings as well as the model arguments (see code for arg options). The fitting script stores the model, the trace, the arguments and the state after adaptation in a date-time stamped folder within the `output` directory. 

The current fitting script does not allow multithreading because the sampling progress cannot be displayed in multithreading mode. Instead, we suggest to open multiple instances of julia to run the script several times (with different seeds). In our experience, running four instances in parallel does not impact performance (but this depends of course on the performance of your machine, please check that the number of instances does not exceed the number of available cores on your machine).  

The following settings/seeds **reproduce** our fitting results:

| model | args | seeds |
| :--- | :--- |:--- |
| Mises 2010-2022 | --betafunc=mises --symb --verbose |25,95,104,222,303,462,1886,1983,9333,68136|
| Mises 2010-2019 | --betafunc=mises --symb --verbose --endyear=2019 --target-acceptance=0.95 |95,104,222,303,462,1886,9333,68136|
| Cosine 2010-2022 | --betafunc=cosine --symb --verbose | 25,644,874,7425,8319,336447|
| Mises 2010-2022, without AR data | --betafunc=mises --symb --verbose --noar |25,95,104,222,303,462,1886,9333,68136,83411 |


<br/> 

To run the **diagnostics/visualization** of the fitting results, run the following from the console:
```sh
julia --project scripts/visualize.jl --verbose --theme=ggplot2 output/foldername
```

where *foldername* is the name of the folder containing the fitting results. You can also combine the results from different fits with the same settings with 
```sh
julia --project scripts/visualize.jl --verbose --theme=ggplot2 --prefix=folderprefix
```
where *folderprefix* is the prefix of the folders with the fits to combine

<br/> 

To generate the **posterior predictive** of the fit, run the following from the console:
```sh
julia --project scripts/postprocessing.jl --verbose output/foldername
```
where *foldername* is the name of the folder containing the fitting results. You can combine the results from different fits as above. 

<br/> 

To simulate the **vaccination scenarios**, run the following from the console:
```sh
julia --project scripts/vaccsim.jl seed output/foldername --verbose 
```
where *seed* is the random number see and *foldername* is the name of the folder containing the fitting results. You can combine the results from different fits as above. 


<br/> 

## Contact
Do you have questions about the study or the code? Contact [Fabienne Krauer](https://www.lshtm.ac.uk/aboutus/people/krauer.fabienne) or [Tor Fjelde](http://www.eng.cam.ac.uk/profiles/tef30) by email or open a GH issue. 

