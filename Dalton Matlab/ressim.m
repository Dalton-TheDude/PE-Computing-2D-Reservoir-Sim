% Formation and cell dimensions
x=50000;            %feet
y=50000;            %feet
thick=20;           %feet
nx=3;              %number of cells on the x axis
ny=3;              %number of cells on the y axis
dx= x ./ (nx);        %feet/cell
dy= y ./ (ny);        %feet/cell

% Time Stuff
ndays=10;           %days
dt=1;               %time step size in days
nsteps=ndays*dt;    %total number of time steps

% Formation properties
p0=7000;            %psi
bhp=2000;           %psi
cr=0.00005;        %1/psi
cf=0.0005;         %1/psi
visc=250;           %cP
k0=200;             %mD
phi0=0.4;           
b0=1.2;             %resb/stb

% Well stuff
wellinfo= [3, 200]; % Location, STB/d          
iter=100;           %number of iterations allowed per time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cell List
ncells = nx .* ny;
cell_list = [1:ncells];

% Faces
nfaces = ((nx - 1) .* ny) + ((ny - 1) .* nx);
face_list = [1:nfaces];

P = zeros(1, ncells);
Pold = zeros(1, ncells);
kinit = zeros(1, ncells);

for i = 1: ny
    for j = 1: nx
        P(1, ((i-1) .* ny) + j) = p0;
        Pold(1, ((i-1) .* ny) + j) = p0 + 100;
        kinit(1, ((i-1) .* ny) + j) = k0;
    end
end


[residx, jacobx] = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, x, y, nx, ny, nfaces)



