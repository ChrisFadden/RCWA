function [ src ] = setupSrc( fileName,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

degrees = pi / 180;

src.theta = 60 * degrees;   %Elevation Incident Angle
src.phi = 30*degrees;       %Azimuthal Incident Angle
src.pte = 1/sqrt(2);                %Source Polarization
src.ptm = 1i/sqrt(2);                %Source Polarization

%Polarization unit vector
n = [0; 0; 1];
k = [sin(src.theta)*cos(src.phi); sin(src.theta)*sin(src.phi); cos(src.theta)];

if(src.theta == 0)
   src.ate = [0; 1; 0]; 
else
    %multiply by -1 to match bench...   check this
    src.ate = -1 * cross(n,k) ./ norm(cross(n,k));
end

src.atm = cross(src.ate,k) / norm(cross(src.ate,k));

%Compute source vector
delta = zeros(N,1);
delta(ceil(N/2)) = 1;

Esrcx = src.ate(1)*src.pte + src.atm(1)*src.ptm;
Esrcy = src.ate(2)*src.pte + src.atm(2)*src.ptm;

src.esrc = [Esrcx .* delta; Esrcy .* delta];


end

