function [Im]=x2k(ik)

% x2k
%
% Data in x-space to k-space
%
% INPUTS:
%	ik:     x-space data
%
% Note: Data format coresponds to method 1 in k2x
%
% PARALLEL MRI TOOLBOX
%
% Santiago Aja-Fernandez, LPI
% www.lpi.tel.uva.es/~santi
% Valladolid, 28/05/2012

if ndims(ik)==2
    Im=fftshift(fftshift(fft2(ik),1),2);
elseif ndims(ik)==3
    [Mx,My,Mz]=size(ik);
    for jj=1:Mz
       Im(:,:,jj)=fftshift(fftshift(fft2(ik(:,:,jj)),1),2); 
    end
elseif ndims(ik)==4
    [Mx,My,Mz, Mu]=size(ik);
    for jj=1:Mz
    for ii=1:Mu
       Im(:,:,jj,ii)=fftshift(fftshift(fft2(ik(:,:,jj,ii)),1),2); 
    end
    end
else
    error('Wrong dim');
end
end

