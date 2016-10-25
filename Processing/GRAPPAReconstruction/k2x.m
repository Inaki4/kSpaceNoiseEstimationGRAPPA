function [Im]=k2x(ik,method)

if method==1
	Im=ifft2(fftshift(fftshift(ik,1),2));
elseif method==2
	Im=ifftshift(ifft2(ifftshift(ik)));
elseif method==3
	Im=fftshift(ifft2(ifftshift(ik)),2); 
end

