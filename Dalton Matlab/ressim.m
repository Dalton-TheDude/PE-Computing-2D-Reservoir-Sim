clear all

% Formation and cell dimensions
nx =        25;              %number of cells on the x axis
ny =        25;              %number of cells on the y axis
dx =        100;      %feet/cell
dy =        100;      %feet/cell
ncells =    nx .* ny;

% Time Stuff
time =      1;
t_final =   10;             %days
dt =        0.01;              %time step size in days
nsteps =    t_final ./ dt;   %total number of time steps
iter =      15;             %number of iterations allowed per time step
tol =       0.01;           % Tolerance in each cell
sum_tol =   tol .* ncells;  % Sum of the tolerances in each cell

% Formation properties
p0 =        3000;           %initial pressure (psi)
cr =        1.0e-5;        %1/psi
cf =        1.0e-4;         %1/psi
visc =      2.5;            %cP
k0 =        0.1 + (100 - 0.1) .* rand(nx .* ny,1);            %mD
phi0 =      0.1 + (0.3 - 0.1) .* rand(nx .* ny,1);            %porosity
b0 =        1.2;            %resb/stb

% Well stuff
well_list = [70, 50];       % Location (cell #), STB/d 
Nwells = size(well_list, 1);
Pwf = zeros(nsteps, Nwells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Well List
well_info = zeros(ncells, 1);
if Nwells > 0
    well_info(well_list(:, 1)) = well_list(:, 2);
end

% Cell List
cell_list = [1:ncells];

% Faces
nfaces = ((nx - 1) .* ny) + ((ny - 1) .* nx);
face_list = [1:nfaces];

P = p0 .* ones(ncells, 1);
Pold = p0 .* ones(ncells, 1);
kinit = k0 .* ones(ncells, 1);


conn_list = connection_list(nx, ny);


[resid, jacob] = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info);

count = 1;
 while (time <= t_final)
     resid = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info);
          
     count = 1;
     while (norm(resid, 2) > tol) && (count < iter)

         [resid, jacob] = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info);
         test2 = jacob\resid;
         norm(resid, 2);
         P = P - (jacob\resid);
         count = count + 1;
     end
     
     for i = 1: Nwells
        Pwf(count, i) = P(well_list(i, 1));
     end
     Pold = P;
     time = time + dt;
     
 end
 
 
     P_plot = zeros(ny, nx);
     count = 1;
     for i = 1: ny
         for j = 1: nx
             P_plot(i, j) = P(count);
             count = count+1;
         end
     end
figure
spy(sparse(jacob))
figure
imagesc(P_plot)
figure
surf(P_plot)
figure
plot(Pwf)
 
 