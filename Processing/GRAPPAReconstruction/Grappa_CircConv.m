function [ Krec ] = Grappa_CircConv( Ksub,N,kernel,R,L)
%GRAPPA_CIRCCONV performs GRAPPA reconstruction as a circular convolution
% Inputs:
%   - Ksub: subsampled k-space image.
%   - N: GRAPPA weights for the interpolation of the missing lines.
%   - kernel: kernel size used for the interpolation.
%   - R: acceleration factor
%   - L: number of coils
% Output:
%   - Krec: reconstructed k-space
%
%	Copyright (C) 2016 Iñaki Rabanillo <irabvil@lpi.tel.uva.es>
%	Laboratorio de Procesado de Imagen, Universidad de Valladolid
%	www.lpi.tel.uva.es
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%First we build the convolution kernels
NW=prod(kernel);
GrappaWeights=zeros([NW*L*(R-1) L]);
N_struct=squeeze(struct2cell(N));
N_struct=N_struct(1,1:(R-1));        
AcqLine=cellfun(@(Cell) strfind(Cell,'x'),N_struct,'UniformOutput',1);
[~,Order]=sort(AcqLine,'descend');
Ind1=repmat((1:(kernel(2)*L))',length(N(1).n)/(kernel(2)*L),1);
for Ri=1:(R-1)
    Ind2=reshape(repmat((kernel(2)*L)*((Ri-1)+(R-1)*(0:length(N(1).n)/(kernel(2)*L)-1)),(kernel(2)*L),1),[],1);
    Ind=Ind1+Ind2;
    GrappaWeights(Ind,:)=N(Order(Ri)).n;
end 

%Now that we have the interpolation weights properly ordered, we can build
%the kernel
if length(size(GrappaWeights))==2
	GrappaWeights=ajusta_coef(GrappaWeights,kernel,R);
end

[X,Y,Z,W]=size(GrappaWeights);

if Z~=W
	error('Wrong nunber of coeffs');
end

Wk_conv=coef_conv(GrappaWeights,kernel,R);

%Now we perform the reconstruction
WkCell=num2cell(Wk_conv,[1 2 3]);
KrecCell=cellfun(@(x) ConvCoil(Ksub,x),WkCell,'UniformOutput',0);
Krec=permute(cell2mat(KrecCell),[1 2 4 3]);

end

%Function to perform the actual coil by coil reconstruction as a
%convolution
function CoilRec=ConvCoil(Ksub,kernelCoil)
    KsubCell=num2cell(Ksub,[1 2]);
    kernelCoilCell=num2cell(kernelCoil,[1 2]);
    KrecCell=cellfun(@(x,y) circconv2d(x,y),KsubCell,kernelCoilCell,'UniformOutput',0); %Why conj(y)? Because conv2 applies conj to the filter h, but since it is a built-in function, it can't be modified, so we correct in advance the error.
    Krec=cell2mat(KrecCell);
    CoilRec=sum(Krec,3);
end

