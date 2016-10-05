function [ device] = setupDevice(fn,grid,NH,device,layer)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    er_fn = 'ER_Layer';
    fext = '.dat';
    er_fn = strcat(fn,er_fn);
    er_fn = strcat(er_fn,num2str(layer));
    er_fn = strcat(er_fn,fext);
    
    ER(:,:,1) = importdata(er_fn);
        
    device.ERC(:,:,layer) = convmat(ER(:,:,1),NH,NH);
    
     ur_fn = 'UR_Layer';
     ur_fn = strcat(fn,ur_fn);
     ur_fn = strcat(ur_fn,num2str(layer));
     ur_fn = strcat(ur_fn,fext);
    
    UR(:,:,1) = importdata(ur_fn);
    device.URC(:,:,layer) = convmat(UR(:,:,1),NH,NH);
end

