# malaria_within_host
- Codes for solving PDEs of a within-between-host malaria model related to the paper ["Immunity can impose a reproduction-survival tradeoff on human malaria parasites"](https://www.biorxiv.org/content/10.1101/2025.01.27.635035v1)
- Codes for analysing fitness of constant and age-of-infection varying transmission investment strategies of parasites

# User guide
- The file **run_full_model.m** runs the model (first the within-host dynamics are resolved and then the human-vector interactions are resolved).
- The parameters are set globally by calling the file **baseline_parameter_set.m**
- Time stepping schemes for the within-host and human-vector dynamics are contained in the functions **within_host_model.m** and **human_vector_model.m** respectively.
- To generate the fitness landscape for fixed (constant) transmission investment (c) strategies run the script **run_withinhost_constant_strats.m**
- For optimizing nonconstant strategies, i.e. transmission investment c is a function of age of infection, use **run_strategy_optimization.m** being careful to ensure that the spline and discretization settings in this script match those in **withinhost_model_optimization.m**