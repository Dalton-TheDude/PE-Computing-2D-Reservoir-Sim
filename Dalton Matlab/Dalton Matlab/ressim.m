clear all

% Formation and cell dimensions
x =         5000;            %feet
y =         5000;            %feet
nx =        5;              %number of cells on the x axis
ny =        5;              %number of cells on the y axis
dx =        x ./ (nx);      %feet/cell
dy =        y ./ (ny);      %feet/cell
ncells =    nx .* ny;

% Time Stuff
time =      1;
t_final =   10;             %days
dt =        1;              %time step size in days
nsteps =    t_final ./ dt;   %total number of time steps
iter =      15;             %number of iterations allowed per time step
tol =       0.01;           % Tolerance in each cell
sum_tol =   tol .* ncells;  % Sum of the tolerances in each cell

% Formation properties
p0 =        4000;           %initial pressure (psi)
cr =        0.00005;        %1/psi
cf =        0.0005;         %1/psi
visc =      250;            %cP
k0 =        100;            %mD
phi0 =      0.4;            %porosity
b0 =        1.2;            %resb/stb

% Well stuff
well_list = [13, 500];       % Location (cell #), STB/d 
Nwells = size(well_list, 1);
Pwf = zeros(nsteps, Nwells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Well List
well_info = zeros(ncells, 1);
if Nwells > 0
    well_info(well_list(:, 1)) = -well_list(:, 2);
end

% Cell List
cell_list = [1:ncells];

% Faces
nfaces = ((nx - 1) .* ny) + ((ny - 1) .* nx);
face_list = [1:nfaces];

P = p0 .* ones(ncells, 1);
Pold = p0 .* ones(ncells, 1);
kinit = k0 .* ones(ncells, 1);

% The loop I've been using to set up the pressures and permeabilities
% for i = 1: ncells
%         P(i, 1) = p0 + (i .* 10);
%         Pold(i, 1) = p0 + 100;
%         kinit(i, 1) = k0;
% end

conn_list = connection_list(nx, ny);


[resid, jacob] = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, x, y, nx, ny, conn_list);

count = 1;
 while (time <= t_final)

      resid = resid + (dt .* well_info);
      b = (P(:, 1) + resid); %(dt .* well_info(:, 1)));
     test1 = jacob\b;
     test2 = jacob\resid;
     count = 1;
     
     while sum(abs(resid(:, 1))) > sum_tol
         if count > iter
            break
         end
         
         for i = 1: Nwells
            Pwf(count, i) = P(well_list(i, 1));
         end
         b = (P(:, 1) + (dt .* well_info(:, 1))) + resid;
         Pold = P;
         P = (jacob\b);
         
         resid = discretize(P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, x, y, nx, ny, conn_list);
         
         count = count + 1;
     end
     time = time + dt;
     
 end
 
figure
spy(sparse(jacob))
 
     P_plot = zeros(ny, nx);
     count = 1;
     for i = 1: ny
         for j = 1: nx
             P_plot(i, j) = P(count);
             count = count+1;
         end
     end
figure
imagesc(P_plot)
figure
surf(P_plot)
figure
plot(Pwf)
 
 