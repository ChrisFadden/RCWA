%*************************************************
% Simulation Parameters
%*************************************************
close all;
clc;
clear all;

%******************
% Source Parameters
%******************
lam0 = 2 * 10^-3;   %Free Space Wavelength
ginc = [0;0;1];     %Incident wave vector
EP = [0;1;0];       %Source Polarization

%***********************
% Trans / Ref Parameters
%***********************
urR = 1.0;  %Permeability of Reflection Region
erR = 2.0;  %Permittivity of Reflection Region
urT = 1.0;  %Permeability of Transmission Region
erT = 9.0;  %Permittivity of Transmission Region

%******************
% Device Parameters
%******************
urd = 1.0;  %Permeability of device
erd = 6.0;  %Permittivity of device

%*****************
% Layer Parameters
%*****************
Lx = 17.5 * 10^(-3);  %period in x-direction
Ly = 15 * 10^(-3);    %period in y-direction
d1 = 5 * 10^(-3);     %thickness of layer 1
d2 = 3 * 10^(-3);     %thickness of layer 2
w = 0.8 * Ly;

%****************
% RCWA Parameters
%****************
Nx = 512;                 %x points in real-space grid
Ny = round(Nx * Ly / Lx); %y points in real-space grid
PQ = 1 * [1, 1];          %number of spatial harmonics for x, y

%***********************************************************************
% Build Device on grid
%***********************************************************************
dx = Lx / Nx;       %grid resolution in x
dy = Ly / Ny;       %grid resolution in y
xa = [0:Nx-1] * dx; %x axis array
xa = xa - mean(xa); %center x array at 0
ya = [0:Ny-1] * dy; %y axis array
ya = ya - mean(ya); %center y array at 0

%Initialize Layers to device
UR = urd * ones(Nx,Ny,2);
ER = erd * ones(Nx,Ny,2);
L = [d1, d2];

%Build Triangle (Layer 1)
h = 0.5 * sqrt(3) * w;
ny = round(h/dy);
ny1 = round((Ny-ny)/2);
ny2 = ny1+ny - 1;

for ny = ny1:ny2
  f = (ny - ny1) / (ny2 - ny1);
  nx = round(f*w/Lx * Nx);
  nx1 = 1 + floor((Nx - nx)/2);
  nx2 = nx1 + nx;
  ER(nx1:nx2,ny,1) = erR;
end

h = surf(ER(:,:,1),'EdgeColor','None','facecolor','interp');
colormap(jet);
colorbar;
view(2)




































