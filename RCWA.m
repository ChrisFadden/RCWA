%*************************************************
% Simulation Parameters
%*************************************************
close all;
clc;
clear all;

%******************
% Source Parameters
%******************
lam0 = 2 * 10^-2;   %Free Space Wavelength
ginc = [0;0;1];     %Incident wave vector
EP = [0;1;0];       %Source Polarization
theta = 0;
phi = 0;
pte = 1;
ptm = 0;


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
Nharmonics = 3;           %number of spatial harmonics for x, y

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Kx, Ky, Kz Matrices <- NOT FINISHED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ninc = sqrt(erR *urR);
kinc = ninc.*[sind(theta)*cosd(phi),sind(theta)*sind(phi),cosd(theta)];
k0 = 2*pi / lam0;

m = 1;
n = 1;
pqm = 1;
pqn = 1;
pq_bound = floor(Nharmonics/2);
pq2 = -pq_bound:pq_bound;

while(m <= Nharmonics^2)
  kx(m) = kinc(1) - pq(pqm)*(2*pi) / (k0*Lx); 
  n = 1;
  pqn = 1;
  while(n <= Nharmonics^2)
    ky(n) = kinc(2) - pq(pqn)*(2*pi) / (k0*Ly);
    kz_ref(m,n) = conj(sqrt(urR*erR - kx(m)^2 - ky(n)^2));
    kz_trn(m,n) = conj(sqrt(urT*erT - kx(m)^2 - ky(n)^2)); 
    if(mod(n,Nharmonics) == 0)
        pqn = pqn+1;
    end
    n=n+1;
  end
  if(mod(pqm,Nharmonics) == 0)
    pqm = 0;
  end
  m = m+1;
  pqm = pqm+1;
end

Kx = diag(kx(:));
Ky = diag(ky(:));
Kz_ref = diag(diag(kz_ref(:,:)));
Kz_trn = diag(diag(kz_trn(:,:)));

%for m = 1:Nharmonics
    %for n = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reflection Side Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i = 1:length(L)
%    ERC = convmat(ER(:,:,i),PQ(1),PQ(2));
%    URC = convmat(UR(:,:,i),PQ(1),PQ(2)); 
%end






































