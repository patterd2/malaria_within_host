function [HS, HI, VS, VI] = human_vector_model(h, t0, X_max, T_max, xV_max, HS0, HI0, VS0, VI0, G)

global P

t = (0:h:T_max)';
nt = length(t);
xV = (0:h:xV_max)';
nxV = length(xV);
x = (t0:h:X_max)';
nx = length(x); % all discretization parameters chosen equal

% allocation of memory
HS = NaN(nt,1); % HS(t) -> HS(mh)
HI = NaN(nt,nx); % HI(t,x) -> HS(mh,nh)
VS = NaN(nt,1); % VS(t)
VI = NaN(nt,nxV); % VI(t,xV)

% assisgn initial conditions
HS(1) = HS0;
HI(1,:) = HI0; % HI(0,x)
VS(1) = VS0;
VI(1,:) = VI0; % VI(0,tauV)

%% time evolution
for m = 1:nt-1 % time stepping
    %fprintf('%i \n',m);
    %keyboard;
    VS(m+1) = (VS(m)/h + RA(m*h))./(1/h + (P.b*h/P.N)*sum( G.*(HI(m,:)') ));
    HS(m+1) = (HS(m)/h)./(1/h + h*P.betaVH*(P.b/P.N)*sum( gammaV(xV,h).*VI(m,:)'));
    HI(m+1,2:end) = HI(m,1:end-1);
    HI(m+1,1) = P.betaVH*(P.b/P.N)*HS(m+1)*h*sum(gammaV(xV,h).*VI(m,:)');
    VI(m+1,2:end) = (VI(m,1:end-1)/h)./(1/h + P.deltaA);
    %keyboard;
    VI(m+1,1) = (P.b/P.N)*VS(m)*h*sum(betaHV(G).*(HI(m,:)'));
end

end
