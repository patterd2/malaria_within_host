# malaria_within_host
Codes for solving PDEs of a within-between-host malaria model

# User guide:
The file run_full_model.m runs the model (first the within-host dynamics are resolved and then the human-vector interactions are resolved).
The parameters are set globally by calling the file baseline_parameter_set.m
Time stepping schemes for the within-host and human-vector dynamics are contained in the functions within_host_model.m and human_vector_model.m respectively.
