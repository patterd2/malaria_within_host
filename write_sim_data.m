% write data to text file
immune_removal = (1-exp(-P.theta*A));
infection_level = h*sum(I,2); % number of iRBCs (integral over tau)
infectiousness = betaHV(G);

mat_data = [G B infection_level A infectiousness immune_removal];
writematrix(x/24,'x_data_days.txt','Delimiter',' ');
writematrix(mat_data,'traj_data_noImmunity.txt','Delimiter',' ');
writestruct(P,'parameter_data_noImmunity.xml');
