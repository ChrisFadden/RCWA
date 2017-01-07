function [Kx, Ky, KzT,KzR] = calcK(src,grid,Nharmonics)

  ninc = sqrt(grid.erR * grid.urR);
  kinc.x = ninc * sin(src.theta)*cos(src.phi);
  kinc.y = ninc * sin(src.theta)*sin(src.phi);
  kinc.z = ninc * cos(src.theta);
  k0 = 2*pi / grid.lam0;
  
  %assume equal and odd number of spatial harmonics in x,y directions
  M = (-(Nharmonics-1)/2):((Nharmonics-1)/2);
  N = length(M)^2; %Size of matrices 
  
  kx = kinc.x * ones(1,sqrt(N));
  ky = kinc.y * ones(1,sqrt(N));
  for m = 1:sqrt(N)
    kx(m) = kinc.x - (2*pi*M(m))/(k0*grid.Lx);
    ky(m) = kinc.y - (2*pi*M(m))/(k0*grid.Ly);
  end

  Kx = zeros(N,N); 
  Ky = zeros(N,N);
  
  for m = 1:sqrt(N)
    for n = 1:sqrt(N)
      ind = (m-1)*sqrt(N) + n;
      Kx(ind,ind) = kx(n);
      Ky(ind,ind) = ky(m);
    end
  end
  
  %CHECK THIS
  KzT = conj(sqrt(grid.erT*grid.urT*eye(N) - Kx.^2 - Ky.^2));
  KzR = -1 * conj(sqrt(grid.erR*grid.urR*eye(N) - Kx.^2 - Ky.^2));

end

