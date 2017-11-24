global params

% Formation and cell dimensions
params.x=50000;            %feet
params.y=50000;            %feet
params.thick=20;           %feet
params.nx=10;              %number of cells on the x axis
params.ny=10;              %number of cells on the y axis
params.dx= x ./ (nx);        %feet/cell
params.dy= y ./ (ny);        %feet/cell

% Time Stuff
params.ndays=10;           %days
params.dt=1;               %time step size in days
params.nsteps=ndays*dt;    %total number of time steps

% Formation properties
params.p0=7000;            %psi
params.bhp=2000;           %psi
params.cr=0.00005;        %1/psi
params.cf=0.0005;         %1/psi
params.visc=250;           %cP
params.k0=200;             %mD
params.phi0=0.4;           
params.b0=1.2;             %resb/stb

% Well stuff
params.wellinfo= [3, 200]; % Location, STB/d          
params.iter=100;           %number of iterations allowed per time step