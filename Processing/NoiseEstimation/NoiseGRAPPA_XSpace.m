function [NoiseMap,NoiseMapXY,NoiseMapCoils] = NoiseGRAPPA_XSpace(GammaK,CK,kernel,L,N,Walsh,Size,Lines,FlagACS,LinesACS)
%NOISEGRAPPA_APROX Generates a noise map based on the approximation of the
%GRAPPA reconstruction as a product in the image space. 
%Input parameters:
%   - GammaK: Covariance matrix for a complex multivariate gaussian (see:
%   https://en.wikipedia.org/wiki/Complex_normal_distribution). 
%   - CK: Covariance matrix for a complex multivariate gaussian (see:
%   https://en.wikipedia.org/wiki/Complex_normal_distribution). 
%   - kernel=[ky kx]: size of the kernel used for the reconstruction.
%   - L: number of coils
%   - N: kernels used for the reconstruction.
%   - Walsh: 3D array [Size,Size,L] with the vector to combine all the
%   coils for each pixel in order to get the final image.
%   - Size of the image.
%   - Lines: acquired lines for the GRAPPA reconstruction
%   - FlagACS: flag indicating if the ACS are present in the acquisition.
%   - LinesACS: acquired lines in the ACS region.
%Output parameters
%   - NoiseMap: complex noise maps for the coil combined image. Contains a
%   complex number for every part whose real/imag part is the variance of
%   the real/imag part for that pixel.
%   - NoiseMapXY: noise maps for the coil combined image containing the 
%   cross correlation between real and imaginary components for every 
%   pixel.
%   - NoiseMapCoils: complex noise maps for every coil previous to 
%   combination. the coil combined image. Contains a complex number for 
%   every part whose real/imag part is the variance of the real/imag part 
%   for that pixel.
%
%	Copyright (C) 2016 Iñaki Rabanillo <irabvil@lpi.tel.uva.es>
%	Laboratorio de Procesado de Imagen, Universidad de Valladolid
%	www.lpi.tel.uva.es
%
%   This is an extension of the noise maps computed from an equivalent
%   image-space reconstruction reported in 
%   "General Formulation for Quantitative G-factor Calculation in GRAPPA 
%   Reconstructions" by F.A. Breuer, S.A.R. Kannengiesser, M. Blaimer, 
%   N. Seiberlich, P.M. Jakob and M. A. Griswold, 
%   Mag. Reson. Med.  62(3):739-746.  July 2009
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

%We obtain the Gamma matrix in the image domain, which is inmediate due to
%the orthogonality of the FFT. The C matrix is more complicated to obatin.
GammaX=GammaK/(Size^2);

%First we identify the regions with different acceleration
NCell=squeeze(struct2cell(N));
PatternCell=NCell(1,:);
RPattern=cell2mat(cellfun(@(x) x(2)-x(1),cellfun(@(x) strfind(x,'*'),PatternCell,'UniformOutput',0),'UniformOutput',0));
R=unique(RPattern);
RLines=diff(Lines);
BoundaryLines=[Lines(1); Lines(1+find(diff(RLines))); Lines(end)];
RRegion=RLines([1;diff(RLines)]~=0);
i=1;
if(FlagACS)
    fi_Ri=length(LinesACS)/Size;
    GammaMapKernel(:,:,:,:,i)=repmat(permute(fi_Ri*GammaX,[3 4 1 2]),[Size Size 1 1]); %If we do it by considering the convolution, the kernel would be a delta centered at 0, whose fourier transform is a constant 1/Size², that would appear twice in the variance analysis, and that would dissapear taking into account thatIt's Size^4 because reconstruction in the image space is X*Y->Size²·(x·y), so variance analysis includes Size² at both sizes X*Y->Size²·(x·y), which again appears twice.
    CMapKernel(:,:,:,:,i) = CXfromK(CK,Size,L,LinesACS); %CX matrix for the ACS lines in the image domain
    i=i+1;
end
for Ri=R
    IndxBoundsRi=find(RRegion==Ri);
    firstLineRi=BoundaryLines(IndxBoundsRi);
    lastLineRi=BoundaryLines(1+IndxBoundsRi);
    numLinesRi=sum(lastLineRi-firstLineRi)/Ri; %+length(IndxBoundsRi);
    if any(((Size/2)<lastLineRi)&((Size/2)>firstLineRi))
        numLinesRi=numLinesRi+1; %If we are in the middle of the k-space, we count both extremes
        LinesRi=firstLineRi:Ri:lastLineRi;
    else
        IndxLow=find(lastLineRi<Size/2);
        IndxUp=find(lastLineRi>=Size/2);
        LinesLow=reshape(cell2mat(arrayfun(@(x,y) x:Ri:(y-Ri),firstLineRi(IndxLow),lastLineRi(IndxLow),'UniformOutput',0)).',[],1);
        LinesUp=reshape(cell2mat(arrayfun(@(x,y) (x+Ri):Ri:y,firstLineRi(IndxUp),lastLineRi(IndxUp),'UniformOutput',0)).',[],1);
        LinesRi=[LinesLow;LinesUp]; %If not, we don't include the boundary closer to the center
    end
    fi_Ri=numLinesRi/Size;
    Ni=N(find(RPattern==Ri));
    [GammaMapKernel(:,:,:,:,i),CMapKernel(:,:,:,:,i)]=NoiseXSpace_Kernel(GammaX,CK,kernel,Ri,L,Ni,Size,LinesRi,fi_Ri);
%     GammaMapKernel(:,:,:,:,i)=fi_Ri*(Size^4)*GammaMapKernelRi; %It's Size^4 because reconstruction in the image space is X*Y->Size²·(x·y), so variance analysis includes Size² at both sides
%     CMapKernel(:,:,:,:,i)=fi_Ri*(Size^4)*CMapKernelRi; %Breuer paper includes Ri multiplying, but this only applies if normalized kernel is used in the equation 
%     NoiseKernel(:,:,i)=Ri*fi*(Size^2)*NoiseXSpace_Kernel(Gamma,C,kernel,Ri,L,Ni,Walsh,Size);
    i=i+1;
end
GammaMap=sum(GammaMapKernel,5);
CMap=sum(CMapKernel,5);
GammaMapCell=permute(num2cell(permute(GammaMap,[3 4 1 2]),[1 2]),[3 4 1 2]);
CMapCell=permute(num2cell(permute(CMap,[3 4 1 2]),[1 2]),[3 4 1 2]);
WalshAux=permute(Walsh,[3 1 2]);
NoiseMapGamma=cellfun(@(x,y) y'*x*y,GammaMapCell,permute(num2cell(WalshAux,[1]),[2 3 1]));
NoiseMapC=cellfun(@(x,y) y'*x*conj(y),CMapCell,permute(num2cell(WalshAux,[1]),[2 3 1]));
NoiseMap=complex(real(NoiseMapGamma+NoiseMapC),real(NoiseMapGamma-NoiseMapC))/2;
NoiseMapXY=imag(-NoiseMapGamma+NoiseMapC)/2;
NoiseMapCoils=reshape(permute(complex(real(GammaMap+CMap),real(GammaMap-CMap))/2,[3 4 1 2]),[L L Size^2]);
IndDiag=cumsum([1:(L+1):L^2; L^2.*ones(Size^2-1,L)])';
NoiseMapCoils=permute(reshape(NoiseMapCoils(IndDiag),[L Size Size]),[2 3 1]);
end

function [GammaMapKernel,CMapKernel]=NoiseXSpace_Kernel(GammaX,CK,kernel,R,L,N,Size,Lines,fi_Ri)
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

    Wx=coef_x(GrappaWeights,R,kernel,Size,1,1);
    WxAux=permute(Wx,[3 4 1 2]);
    GammaMapKernel=fi_Ri*(Size^4)*permute(cell2mat(permute(cellfun(@(x) x.'*GammaX*conj(x),permute(num2cell(WxAux,[1 2]),[3 4 1 2]),'UniformOutput',0),[3 4 1 2])),[3 4 1 2]);
    CX = CXfromK(CK,Size,L,Lines); %CX matrix for the ACS lines in the image domain
    CXCell=permute(num2cell(permute(CX,[3 4 1 2]),[1 2]),[3 4 1 2]);
    CMapKernel=(Size^4)*permute(cell2mat(permute(cellfun(@(x,y) x.'*y*x,permute(num2cell(WxAux,[1 2]),[3 4 1 2]),CXCell,'UniformOutput',0),[3 4 1 2])),[3 4 1 2]);
%     NoiseMapGamma=cellfun(@(x,y) y'*(x'*Gamma*x)*y,permute(num2cell(WxAux,[1 2]),[3 4 1 2]),permute(num2cell(WalshAux,[1]),[2 3 1]));
%     NoiseMapC=cellfun(@(x,y) y'*(x'*Gamma*conj(x))*conj(y),permute(num2cell(WxAux,[1 2]),[3 4 1 2]),permute(num2cell(WalshAux,[1]),[2 3 1]));
%     NoiseMap=complex(real(NoiseMapGamma+NoiseMapC),real(NoiseMapGamma-NoiseMapC))/2;
end