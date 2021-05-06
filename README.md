# Dengue vaccine trial survival model

This repository contains all code for the Bayesian survival analysis of the efficacy profile of 
the Sanofi-Pasteur CYD-TDV dengue vaccine
(Dengvaxia), developed primarily by Daniel J. Laydon and Neil M. Ferguson of the MRC Centre
for Global Infectious Disease Analysis at Imperial College, London.


## IMPORTANT NOTES

:warning: We are not at liberty to share the underlying individual-level data from each trial, as 
it is proprietary and belongs to Sanofi-Pasteur. 
However, we have provided simulated trial data that approximately preserves case-counts
across multiple strata. When our model is run on this simulated
data, our analysis is largely reproduced. 

:warning: Parameter files and code are named quite esoterically and may be difficult to interpret. 
Updating the code to make it more easily interpretable and user-friendly is a work-in-progress.

## Building

The model is written in C++, using OpenMP to improve performance on multi-core processors. 
C++ source and header files are located in the subdirectory "DengVaxSurvival", 
and R output processing code can be found in the subdirectory "R". 


## Simulated Sample Data

The directory [DengVaxSurvival/Data](./DengVaxSurvival/Data) contains simulated sample data 
for our default model, and for models that separate cases by disease severity, as well as models 
without serotype effects. 

## Running

The target executable of the C++ source code is named `D_MCMC.exe`, which takes a single 
command line argument specifying a parameter file. So the syntax for running the executable is 
`D_MCMC.exe Params_ParticularModelSettings.txt`. 
The code is heavily parallelized and (broadly) will run faster with more cores, ideally on a 
high-performance cluster. However the model can also be run locally if the source code is 
recompiled having turned off/commented out two lines in the header file Macros.h, namely 
`#define USE_CLUSTER` and `#define USE_COMMAND_LINE`. 
If run locally, have set number of cores at 6 (in main.cpp), which again can be changed if desired. 

Note to run on a cluster, you will need to have the .dll files included in [dll_files](./dll_files) copied to whichever directory the
executable (.exe) is run from. 

## Parameter files

The executable reads in a parameter file, named "Params_ParticularModelSettings.txt", 
that govern which features are turned on or off,
(initial) model parameter values, and the outputs tha are required,
(e.g. MCMC chains, augmented data, survival curves, attack rates, hazard ratios, seroprevalence etc.).
Example parameter files are provided in the directory [ParamFiles](./ParamFiles).

By default, the features specified in the parameter files 
produce an output string that identifies the features present in that model run. 
This output string is then included within all output file names. 
The output strings are long, unweildly and difficult to interpret. We apologise for the inconvenience
and are working to make this less painful. 

Our main model (with both age and serotype effects) has output string
"VAC_SILENT_PASSIVE_nENWX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU". 
Model with age effects only has output string 
"VAC_SILENT_PASSIVE_nENWX_FSKs3_AGSVEhetero_AS_Hazmult_fAdjHaz_prs1_2_SFU"
Model with serotype effects only has output string 
"VAC_SILENT_PASSIVE_nENWX_SS_VEs_FSKs3_fAdjHaz_prs1_2_SFU"
Model with neither age nor serotype effects has output string 
"VAC_SILENT_PASSIVE_nENWX_FSKs3_fAdjHaz_prs1_2_SFU".

The above models all assume the "vaccine-as-silent-infection" mechanism, and do not separate cases by disease severity. 
The parameter file, 
"Params_K_SEROPOS_PASSIVE_nENWX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU.txt", 
models a variant without such immune priming. 

Finally, parameters that model hospitalised and non-hospitalised disease separately are given in 
file "Params_VAC_SILENT_PASSIVE_nENW_MILDSEVERE_hospX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU.txt".  
Modelling severe and non-severe disease separately can be accomplished using 
"Params_VAC_SILENT_PASSIVE_nENW_MILDSEVEREX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU.txt". 
Please note that these parameter files use different simulated data files "SimData_Hosp.txt" and "SimData_Severe.txt"
Default output strings are created in the C++ function 
StructureDefs.h::Housekeeping_Struct::CreateOutputString, and in the R function DirectoriesEtc.R::ChooseOutputString

## Output

All model outputs are `.txt` files. If run locally, model output (e.g. parameter chains, survival tables, attack rates, estimated 
age-specific seroprevalence) is by default stored in 
[DengVaxSurvival/Output](./DengVaxSurvival/Output), (although all outputs are untracked by git) as 
they can be quite large. 

The folder [R](./R) contains scripts to process and plot model output. By default, these scripts
expect model output to be located in [DengVaxSurvival/Output](./DengVaxSurvival/Output).


### Relevant papers

The following papers are relevant to the model and trial data. Please note that some of them
may require a subscription.

- <https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(12)61428-7/fulltext>
- <https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(14)61060-6/fulltext>
- <https://www.nejm.org/doi/10.1056/NEJMoa1411037?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0www.ncbi.nlm.nih.gov>
- <https://www.nejm.org/doi/full/10.1056/nejmoa1506223>
- <https://science.sciencemag.org/content/353/6303/1033.abstract>

## Copyright and Licensing

Source code is licensed under the GPLv3.

It is Copyright Imperial College of Science, Technology and Medicine. The
lead developers are Daniel J. Laydon and Neil M. Ferguson. 

This repository includes code modified from
[RANLIB](https://people.sc.fsu.edu/~jburkardt/c_src/ranlib/ranlib.html) which
is licensed under the LGPLv3.

