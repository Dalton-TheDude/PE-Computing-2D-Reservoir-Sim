%Still need:
%flow across each face
%V of fluids in res
%Pressure update
%qcell=sum(flow across each face)
%qcell depends on the distance to the nearest well
%volume of fluid in each cell/whole reservoir plot

x=50000;            %feet
y=50000;            %feet
thick=20;           %feet
nx=50;              %number of cells on the x axis
ny=50;              %number of cells on the y axis
dx=x/(nx-1);        %feet/cell
dy=y/(ny-1);        %feet/cell
ndays=10;           %days
dt=1;             %time step size in days
nsteps=ndays*dt;    %total number of time steps
p0=7000;            %psi
bhp=2000;           %psi
cr=0.000005;        %1/psi
cf=0.00005;         %1/psi
visc=250;           %cP
k0=200;             %mD
phi0=0.4;           
b0=1.2;             %resb/stb
nwells=1;           
iter=100;           %number of iterations allowed per time step
vres=x*y*thick;     %feet^3
vcell=dx*dy*thick;  %feet^3 per cell

p=zeros(nx,ny);     
k=zeros(nx,ny);     
phi=zeros(nx,ny);   
b=zeros(nx,ny);     
y=zeros(nx,ny);     
z=[nx,ny,p0];
peq=zeros((nx*ny),1);
qper=zeros(nwells);
qtot=sum(qper);
wellloc=zeros(nwells);
vfluid=zeros(nx,ny);

%Boundary Pressure Conditions
pwest=p0;
peast=p0;
pnorth=p0;
psouth=p0;

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

for i=1:nx
    for j=1:ny
        vfluid(i,j)=(5.615*dx*dy*thick*phi(i,j))/(b(i,j));
    end
end

syms p1
phi1=phi0*exp(cr*(p1-p0));
phip(p1)=diff(phi1,p1);
y1=cf*(p1-p0);
b1=b0/(1+y1+(0.5*(y1^2)));
bp=diff(b1,p1);

syms pi pii pj pjj
xtemp1=k(i,j)/(visc*b(i,j));
xtemp2=(pii-pi)/2;
xterm=xtemp1*xtemp2/2;

ytemp1=k(i,j)/(visc*b(i,j));
ytemp2=(pjj-pj)/2;
yterm=ytemp1*ytemp2/2;

positionterm=xterm+yterm;
diff(xtemp2, pi)



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

