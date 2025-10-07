# malaria_within_host
- Codes for solving PDEs of a within-between-host malaria model related to the paper ["Immunity can impose a reproduction-survival tradeoff on human malaria parasites"](https://www.biorxiv.org/content/10.1101/2025.01.27.635035v1)
- Codes for analysing fitness of constant and age-of-infection varying transmission investment strategies of parasites

# User guide
- The file **run_full_model.m** runs a single simulation of the within-host model.
- The parameters are set globally by calling the file **baseline_parameter_set.m**
- Time stepping schemes for the within-host dynamics are contained in the function **within_host_model.m**.
- To generate the fitness landscape for fixed (constant) transmission investment (c) strategies run the script **run_withinhost_constant_strats.m**
- For optimizing nonconstant/age-varying strategies, i.e. transmission investment c is a function of age of infection, use **run_strategy_optimization.m** being careful to ensure that the spline and discretization settings in this script match those in **withinhost_model_optimization.m**
- Splines are imported for various mesh spacings from text files labelled 'basisMatrix...txt' and new sets of splines can be generated using **CodeToGenerateBasisMatrix.R**
- The folder 'sensitivity_results' contains optimal strategy weights for a range of different parameter regimes.