function [ device,g] = setupDevice( fileName,g,Nx,Ny,Nharmonics)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%Reflection Region
g.urR = 1.0;          %Permeability of Reflection Region
g.erR = 2.0;          %Permittivity of Reflection Region

%Transmission Region
g.urT = 1.0;          %Permeability of Transmission Region
g.erT = 9.0;          %Permittivity of Transmission Region

%Device region
urd = 1.0;          %Permeability of device
erd = 6.0;          %Permittivity of device

%Initialize Layers to device
UR = urd * ones(Nx,Ny,2);
ER = erd * ones(Nx,Ny,2);

%Build Device (Layer 1)
w = 0.8 * g.Ly;
h = 0.5 * sqrt(3) * w;
ny = round(h/g.dy);
ny1 = round((Ny-ny)/2);
ny2 = ny1+ny - 1;

Lx = 17.5 * 10^(-3);

for ny = ny1:ny2
  f = (ny - ny1) / (ny2 - ny1);
  nx = round(f*w/Lx * Nx);
  nx1 = 1 + floor((Nx - nx)/2);
  nx2 = nx1 + nx;
  ER(nx1:nx2,ny,1) = g.erR;
end

for layer = 1:length(g.L)
    device.ERC(:,:,layer) = convmat(ER(:,:,layer),Nharmonics,Nharmonics);
    device.URC(:,:,layer) = convmat(UR(:,:,layer),Nharmonics,Nharmonics);
end

end

