function [w2]=ajusta_coef(w,kernel,R)

% Adjust size of matrix of GRAPPA coefficients
% NORMALIZATION
% From 8x48 or 48x8 GRAPPA coefs
% Build [2x3x8x8] matrix
% 
% PARALLEL MRI TOOLBOX
%
% Santiago Aja-Fernandez, LPI
% www.lpi.tel.uva.es/~santi
% Valladolid, 28/05/2012

[M,N]=size(w);
if M>N
w=w.';
end
[M,N]=size(w);
%M=8
%N=8x3x2

w2=zeros([(R-1)*kernel(1),kernel(2),M,M]);

for kk=1:M

  Tw=reshape(w(kk,:),[kernel(2)*M,(R-1)*kernel(1)]);

  w2(:,:,:,kk)=reshape(Tw.',[(R-1)*kernel(1),kernel(2),M]);

end




