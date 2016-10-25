function  [MapW]=sensitivy_map(Sz,Ncoil)

% SENSITIVITY MAP creates aritficial sensitivity maps
%
%   The map is created such as the Sum of Squares is 1
%
%                      sum(MapW.^2,3)=ones(Sz)
%
%  INPUTS:
%               Sz=[Mx,My]   Size of the map
%               Ncoil        Number of coils   (even number)
%
%  OUTPUT
%               MapW     [Mx,My,Ncoils] Sensitivity Map
%
% EXAMPLE
%              MapW=sensitivy_map([256,256],4);
%
% Santiago Aja-Fernandez
% Parallel Imaging Toolbox
% Feb. 2012
% www.lpi.tel.uva.es/~santi

if numel(Sz)~=2
   error('Input Size must be 2D');
end
 if rem(Ncoil,2)==1
   error('Number of coils must be even');
end 

Mx=Sz(1);
My=Sz(2);
ejex=1:Mx;
ejey=1:My;
vX=repmat(ejex',[1,My]);
vY=repmat(ejey,[Mx,1]);
vX=pi/2.*vX./max(vX(:));
vY=pi/2.*vY./max(vY(:));


Theta=0:(2*pi/Ncoil):(2*pi-(2*pi/Ncoil));
Theta=Theta(1:end-1);
for ii=1:Ncoil./2
  if (Theta(ii)<=pi/2)
      N1=vX.*cos(Theta(ii))+vY.*sin(Theta(ii));
      N1=N1./max(N1(:)).*pi/2; 
      MapW(:,:,ii)=cos(N1);
      MapW(:,:,Ncoil./2+ii)=sin(N1);
  else %>pi/2
      N1=(vX).*abs(cos(Theta(ii)))+(pi/2-vY).*sin(Theta(ii));
      N1=N1./max(N1(:)).*pi/2; 
      MapW(:,:,ii)=sin(N1);
      MapW(:,:,Ncoil./2+ii)=cos(N1);
      
  end

end

MapW=MapW./sqrt(Ncoil/2);	



