function [B, M, I, IG, G, A] = within_host_model(h, t0, tfinal, B0, M0, I0, IG0, G0, A0)

global P

x = (t0:h:tfinal)';
nx = length(x); % all discretization parameters chosen equal

% allocation of memory
B = NaN(1,nx); % B(x)
M = NaN(1,nx); % M(x)
I = NaN(nx,nx); % I(x,tau)
IG = NaN(nx,nx); % IG(x,tau)
G = NaN(1,nx); % G(x)
A = NaN(1,nx); % A(x)

% assisgn initial conditions
B(1) = B0;
M(1) = M0;
I(1,:) = I0;
IG(1,:) = IG0;
G(1) = G0;
A(1) = A0;
I(1,1) = (1-P.c)*P.p*B0*M0;
IG(1,1) = P.c*P.p*B0*M0;
%% time since infection evolution
for n = 1:nx-1 % evolving on the time since infection time scale
    % evolve ODEs for red blood cells and merozoites
    %fprintf('%i \n',n);
    B(n+1) = (B(n)/h + P.lambda)./(1/h + P.lambda/P.K + P.p*M(n) + P.mu);
    M(n+1) = (M(n)/h + h*P.beta*sum(gamma_fun(h*x).*I(n,:)'))./( 1/h + P.muM +P.p*B(n));
    % BC for I
    I(n+1,1) = (1-P.c)*P.p*B(n+1)*M(n+1);
    % evolve I
    I(n+1,2:end) = (I(n,1:end-1)/h)./(1/h + P.mu + gamma_fun(h*x(2:end)') + P.sigma.*(1-exp(-P.theta*A(n))) );
    % BC for IG
    IG(n+1,1) = P.c*P.p*B(n+1)*M(n+1);
    % evolve IG
    IG(n+1,2:end) = ((1/h)*IG(n,1:end-1))./(1/h + P.muG + gamma_G(h*x(2:end))' );
    % evolve ODEs for gametocytes and the immune activation level
    G(n+1) = ((1/h)*G(n) + h*sum(gamma_G(h*x).*IG(n,:)'))./(1/h + P.muG);
    A(n+1) = A(n) + h*(phi( h*sum(I(n,:)), P.IT, P.s));
end
end
