clear all

% Formation and cell dimensions
nx =        25;              %number of cells on the x axis
ny =        25;              %number of cells on the y axis
dx =        100;      %feet/cell
dy =        100;      %feet/cell
ncells =    nx .* ny;

% Time Stuff
time =      1;
t_final =   100;             %days
dt =        0.1;              %time step size in days
nsteps =    t_final ./ dt;   %total number of time steps
iter =      15;             %number of iterations allowed per time step
tol =       0.1;           % Tolerance in each cell
sum_tol =   tol .* ncells;  % Sum of the tolerances in each cell

% Formation properties
p0 =        3000;           %initial pressure (psi)
cr =        1.0e-5;        %1/psi
cf =        1.0e-4;         %1/psi
visc =      2.5;            %cP
k0 =        50 .* ones(nx .* ny,1);            %mD
phi0 =      0.3 .* ones(nx .* ny,1);            %porosity
b0 =        1.2;            %resb/stb

% Wells located in cells 183 and 198 in a 20 by 20 reservoir
% This is the way I created the fracture
% for i = 183:1:198
%     k0(i, 1) = 245;
% end



% Well stuff
% This one is the one I used for the fracture on a 20x20 grid
% well_list = [183, 100; 198, 100];       % Location (cell #), STB/d

% This is the one that puts wells in the corners
well_list = [nx .* 2 + 3, 50; nx .* (ny-2) - 2, -50];

% This is the one that produces the checkerboard pattern on a 25x24 grid
% well_list = [(nx+3), -50; (2 .* nx -2), -50; ((2 .* nx) - ((nx+1) ./ 2) +1), -50; (12 .* nx +3), -50; (13 .* nx -2), -50; ((13 .* nx) - ((nx+1) ./ 2) +1), -50; (22 .* nx +3), -50; (23 .* nx -2), -50; ((23 .* nx) - ((nx+1) ./ 2) +1), -50; (6 .* nx) + 8, 100; (6 .* nx) + 18, 100; (17 .* nx) + 8, 100; (17 .* nx) + 18, 100]; 

Nwells = size(well_list, 1);

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


P = p0 .* ones(ncells, 1);
Pold = p0 .* ones(ncells, 1);
kinit = k0;


conn_list = connection_list(nx, ny);


count = 1;
 while (time <= t_final)
     resid = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info);
          
     count = 1;
     while (norm(resid, 2) > tol) && (count < iter)

         [resid, jacob] = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info);
         P = P - (jacob\resid);
         count = count + 1;
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