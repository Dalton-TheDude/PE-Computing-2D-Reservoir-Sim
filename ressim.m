x=50000;            %feet
y=50000;            %feet
nx=60;              
ny=50;              
dx=x/(nx-1);        %feet/section
dy=y/(ny-1);        %feet/section
ndays=10;           %days
dt=0.1;             %time step size
nsteps=ndays*dt;    %total number of time steps
p0=5000;            %psi
cr=0.000005;        %1/psi
cf=0.00005;         %1/psi
visc=250;           %cP
k0=200;             %mD
phi0=0.4;           
b0=1.2;             %resb/stb
nwells=1;           
iter=100;           %number of iterations allowed per time step

p=zeros(nx,ny);     
k=zeros(nx,ny);     
phi=zeros(nx,ny);   
b=zeros(nx,ny);     
y=zeros(nx,ny);     
z=[nx,ny,p0];       

%Boundary Pressure Conditions
pwest=2000;
peast=2000;
pnorth=2000;
psouth=2000;

for j=1:ny
    p(1,j)=pnorth;
    p(nx,j)=psouth;
end

for i=1:nx
    p(i,1)=pwest;
    p(i,ny)=peast;
end


%Initial Pressure Conditions
for i=2:(nx-1)
    for j=2:(ny-1)
        p(i,j)=p0;
    end
end

%Permeability in each cell in mD
for i=1:nx
    for j=1:ny
        k(i,j)=k0;
    end
end

%Initial Porosity in each cell
for i=1:nx
    for j=1:ny
        phi(i,j)=phi0;
    end
end

%Initial Formation Volume Factor in each cell
for i=1:nx
    for j=1:ny
        b(i,j)=b0;
    end
end




%Simulation part
%Time
for t=0:dt:ndays
    
    %Updating pressure
    for i=2:(nx-1)
        for j=2:(ny-1)
            p(i,j)=p(i,j)-10;
        end

        %Updating porosity
        for i=1:nx
            for j=1:ny
                phi(i,j)=phi0*exp(cr*(p(i,j)-p0));
            end
        end

        %Updating y
        for i=1:nx
            for j=1:ny
                y(i,j)=cf*(p(i,j)-p0);
            end
        end

        %Updating formation volume factor
        for i=1:nx
            for j=1:ny
                b(i,j)=b0/(1+y(i,j)+(0.5*(y(i,j)^2)));
            end
        end
    end
    surf(p)
    drawnow
end

