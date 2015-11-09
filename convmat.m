function C = convmat(A,P,Q,R)
  %This creates the convolution matrices,
  %used in RCWA and PWEM methods.
  
  [Nx,Ny,Nz] = size(A);
  
  %Default Arguments
  if (nargin == 2)
    Q = 1;
    R = 1;
  elseif (nargin == 3)
    R = 1;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute Number of Harmonics
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Nharmonics = P * Q * R;

  p = [-floor(P/2):+floor(P/2)];
  q = [-floor(Q/2):+floor(Q/2)];
  r = [-floor(R/2):+floor(R/2)];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Find Fourier Coefficients of A
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  A = fftshift(fftn(A)) / (Nx*Ny*Nz);
  
  %%%%%%%%%%%%%%%%
  % Array Indicies
  %%%%%%%%%%%%%%%%
  p0 = floor(Nx/2) + 1; %Added 1
  q0 = floor(Ny/2) + 1; %because of 
  r0 = floor(Nz/2) + 1; %MATLAB indexing...

  %%%%%%%%%%%%%%%%%%%%%
  % Loop Over Harmonics
  %%%%%%%%%%%%%%%%%%%%%
  for rrow = 1:R
    for qrow = 1:Q
      for prow = 1:P
        row = ((rrow-1)*Q*P) + (qrow-1)*P + prow;
        for rcol = 1:R
          for qcol = 1:Q
            for pcol = 1:P
              col = (rcol - 1)*Q*P + (qcol - 1)*P + pcol;
              pfft = p(prow) - p(pcol);
              qfft = q(qrow) - q(qcol);
              rfft = r(rrow) - r(rcol);
              C(row,col) = A(p0+pfft, q0+qfft, r0+rfft); 
            end %pcol
          end %qcol
        end %rcol
      end %prow
    end %qrow
  end %rrow
end %function







