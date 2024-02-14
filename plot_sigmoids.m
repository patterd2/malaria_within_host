global P
grid = 0:0.1:10;
P.IT = 2; 
P.s = 0.2;
phi_vals = phi(grid,P.IT,P.s);
figure(8);
plot(grid,phi_vals,'r','LineWidth',3);
hold on;
xline(P.IT,'--k','LineWidth',3);
%xticks([0 10^3 10^6 10^9]);
axis tight;