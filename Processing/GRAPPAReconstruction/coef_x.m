function Wx=coef_x(coef,tasa,kernel,N,mode,full)

% GRAPPA convolution matrix in the x-domain
%
% INPUTS:
%
% 	coef: GRAPPA coefs
% 	tasa: Acceleration Rate
% 	N: Number of points (NxN size of image)
% 	mode: Reconstruction mode (see k2x)
% 	full: full iDFT
%
%
% PARALLEL MRI TOOLBOX
%
% Santiago Aja-Fernandez, LPI
% www.lpi.tel.uva.es/~santi
% Valladolid, 28/05/2012
%------------------------------------------
%-
if length(size(coef))==2
	coef=ajusta_coef(coef,kernel,tasa);
end

[X,Y,Z,W]=size(coef);

if Z~=W
	error('Wrong nunber of coeffs');
end

Wk_conv=coef_conv(coef,kernel,tasa);

Wx=zeros(N,N,Z,W);
    
%The kernel should be centered in the image
WkFull=zeros(N,N,Z,W);
[hx,hy,~,~]=size(Wk_conv);
cx=ceil(hx/2);cy=ceil(hy/2);
WkFull((N/2-cx+2):(N/2+hx-cx+1),(N/2-cy+2):(N/2+hy-cy+1),:,:)=Wk_conv;

%Now we transform to x-space
Wx=k2x(WkFull,mode);
end