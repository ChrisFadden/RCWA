clear;
clc;

grid = importdata('Grid.dat');
urd = 1.0;
erd = 6.0;
erd2 = 2.0;

%Initialize Layers to device
UR = urd * ones(grid(4),grid(5));
ER = erd * ones(grid(4),grid(5));

%Build Device (Layer 1)
w = 0.8 * grid(3);
h = 0.5 * sqrt(3) * w;

nyFac = (grid(3) / grid(5));

ny = round(h/nyFac);
ny1 = round((grid(5)-ny)/2);
ny2 = ny1+ny - 1;

Lx = 17.5 * 10^(-3);

save ER_Layer2.dat ER;

for ny = ny1:ny2
  f = (ny - ny1) / (ny2 - ny1);
  nx = round(f*w/grid(2) * grid(4));
  nx1 = 1 + floor((grid(4) - nx)/2);
  nx2 = nx1 + nx;
  ER(nx1:nx2,ny) = erd2;
end

save ER_Layer1.dat ER;

save UR_Layer1.dat UR;
save UR_Layer2.dat UR;




