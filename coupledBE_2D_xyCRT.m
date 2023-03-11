%% This program is for solving the 2D coupled burgers' equation
%% A Matlab code written by and developed by HUSSEIN A. H. Muhammed. March 2023.
%% B.Sc. and M.Sc. (Honuors) degree in geophysics, 2010-2022.

%coupled BE : 
% du/dt + u(du/dx) + v(du/dy) = v[(d^2u/dx^2)+(d^2u/dy^2)]
% dv/dt + u(dv/dx) + v(dv/dy) = v[(d^2u/dx^2)+(d^2u/dy^2)]

clear;

%% Set equation variables and parameters
nx = 41;
ny = 41;
nt = 120;
c = 0.1;
dx = 2/(nx-1);
dy = 2/(ny-1);
sigma = 0.009;
nu = 0.001;
dt = sigma*dx*dy/nu;


%% create the computational domain
x = linspace(0,2,nx);
y = linspace(0,2,ny);

%% create the viscosity and velocity matrices
u= ones(ny,nx);
v= ones(ny,nx);
un=u;
vn=v;
comb=ones(ny,nx);

%% intial conditions

u(0.5/dy:1/dy+1, 0.5/dx:1/dx+1) = 2;
u(0.5/dy:1/dy+1, 0.5/dx:1/dx+1) = 2;


[X,Y] = meshgrid(x,y);
% surf(x,y,u)    % to test if the mesh metrices are generated

%% Finite-differnce spatial & temporal discritization loop

for n=1:nt+1
    un = u;
    vn = v;
    for i = 2:(ny-1)
        for j = 2:(nx-1)
        u(1,j) = u(i,j) - (dt/dx) * u(i,j)*(u(i,j) - u(i-1,j)) - (dt/dy) * v(i,j) * (u(i,j)-u(i,j-1)) + (nu+dt/dx^2) * (u(i-1,j)-2*u(i,j)...
            +u(i-1,j)+(nu*dt/dy^2)*u(i,j+1)-2*u(i,j)+u(i,j-1));
        v(1,j) = v(i,j) - (dt/dx) * u(i,j)*(v(i,j) - v(i-1,j)) - (dt/dy) * v(i,j) * (v(i,j)-v(i,j-1)) + (nu+dt/dx^2) * (v(i-1,j)-2*v(i,j)...
            +v(i-1,j)+(nu*dt/dy^2)*v(i,j+1)-2*v(i,j)+v(i,j-1));

        %boundary conditions

        u(1:ny,1)=1;
        u(1,1:nx)=1;
        u(1:nx,ny)=1;
        u(ny, 1:nx)=1;
        v(1:ny,1)=1;
        v(1,1:nx)=1;
        v(1:nx,ny)=1;
        v(ny, 1:nx)=1;
        end

    end
    surf(x,y,u);
    colorbar; 
    c = colorbar;
    c.Label.String = 'Velocity Field (unit/sec.)';
    set(gca, 'ydir', 'norm');
    title('Coupled2DBurgerEquation');
    pause(0.1);


end

