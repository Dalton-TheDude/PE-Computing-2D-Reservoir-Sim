clear all
clc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This solves the diffusivity equation in 2D
% Dp/Dt - K D2P/DX2 = Q
% q(x,0) = q(x,L) = q(0,y) = q(L,y) =0
%     { Q_inj     if (x,y) = (0,0)
% Q = { Q_prod    if (x,y) = (0,0)
%     { 0         otherwise
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT PARAMETERS
K      = 1;   % hydraulic diffusivity
L      = 1;   % length of square domain
Q_inj  = 0;  % injector rate
Q_prod = 10;  % producer rate
PINIT  = 4.0;   % Initial reservoir pressure

% SIMULATION PARAMETERS
N      = 100; % Number of simulation gridblocks
DT     = 0.01; % time-step size
Tfinal = 0.5; % final simulation time to stop at

% SETUP
Nvar = N * N;
h    = L/N;
P    = PINIT * ones( Nvar, 1 );
t    = 0; 

% build coefficient matrix A outside loop since it doesnt change
c = K*DT/h/h;
l = ones(Nvar,1);
lft = l;
lft(N:N:end) = 0;
rgt = l;
rgt(N+1:N:end) = 0;
cntr = -4*l;
cntr(1:N) = cntr(1:N) + 1;
cntr((N-1)*N+1:end) = cntr((N-1)*N+1:end) + 1;
cntr(N:N:end) = cntr(N:N:end) + 1;
cntr(1:N:end) = cntr(1:N:end) + 1;
A = speye(Nvar, Nvar) - ...
    c .* spdiags( [ l, lft, cntr, rgt, l], [-N,-1,0,1,N], Nvar, Nvar );
figure
spy(A)

% TIMESTEPPING LOOP
figure
iter = 0;
while ( t < Tfinal )
    % build A and b
    b = P + DT * [ Q_inj; zeros(Nvar-2,1); -Q_prod];
    
    % solve for new pressure
    P = A\b;
    t = t + DT;
    iter = iter + 1
    
    % PLOT 1D Pressure profile along diagonal between wells
    plot( P(1:N+1:end) )
    hold on

end

% PLOT FINAL ANSWER IN 2D
figure
imagesc( reshape(P,N,N) )
figure
contour( reshape(P,N,N), 1000 )



