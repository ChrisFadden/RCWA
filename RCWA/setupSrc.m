function [ src ] = setupSrc( fileName,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

degrees = pi / 180;

src.theta = 0 * degrees;   %Elevation Incident Angle
src.phi = 0 * degrees;       %Azimuthal Incident Angle
src.pte = 1;               %Source Polarization
src.ptm = 0;               %Source Polarization

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
esrcShape = zeros(sqrt(N),sqrt(N));
center = ceil(sqrt(N)/2);

%Plane Wave Example
esrcShape(center,center) = 1;

%Sine Wave Example
%esrcShape(center-1,center-1) = 1;
%esrcShape(center+1,center+1) = 1;

esrcShape = reshape(esrcShape,[N,1]);

%normalize
esrcShape = esrcShape / norm(esrcShape);

Esrcx = src.ate(1)*src.pte + src.atm(1)*src.ptm;
Esrcy = src.ate(2)*src.pte + src.atm(2)*src.ptm;

src.esrc = [Esrcx .* esrcShape; Esrcy .* esrcShape];
end

