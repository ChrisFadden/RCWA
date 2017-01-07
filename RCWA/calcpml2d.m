function [ sx,sy ] = calcpml2d(NGRID,NPML)
    
    %Constants
    a_max = 3;
    sig_max = 1;
    p = 3;
    eta = 376.73;
    
    sx = ones(NGRID(1),NGRID(2));
    sy = ones(NGRID(1),NGRID(2));
    
    %x-low boundary
    a = zeros(1,NPML(1));
    sig = zeros(1,NPML(1));
    for i = 1:NPML(1)
        a(i) = 1 + a_max * (i/NPML(1))^p;
        sig(i) = sig_max * (sin(pi*i / (2*NPML(1))))^2;
        sx(NPML(1)-i+1,:) = a(i) * (1 + 1j*eta*sig(i));
    end
    
    %x-high boundary
    a = zeros(1,NPML(2));
    sig = zeros(1,NPML(2));
    for i = 1:NPML(2)
        a(i) = 1 + a_max * (i/NPML(2))^p;
        sig(i) = sig_max * (sin(pi*i / (2*NPML(2))))^2;
        sx(NGRID(1) - NPML(2) + i,:) = a(i) * (1 + 1j*eta*sig(i));
    end
    
    %y-low boundary
    a = zeros(1,NPML(3));
    sig = zeros(1,NPML(3));
    for i = 1:NPML(3)
        a(i) = 1 + a_max * (i/NPML(3))^p;
        sig(i) = sig_max * (sin(pi*i / (2*NPML(3))))^2;
        sy(:,NPML(3)-i+1) = a(i) * (1 + 1j*eta*sig(i));
    end
    
    %y-high boundary
    a = zeros(1,NPML(4));
    sig = zeros(1,NPML(4));
    for i = 1:NPML(4)
        a(i) = 1 + a_max * (i/NPML(4))^p;
        sig(i) = sig_max * (sin(pi*i / (2*NPML(4))))^2;
        sy(:,NGRID(2) - NPML(4) + i) = a(i) * (1 + 1j*eta*sig(i));
    end
    
    
end

