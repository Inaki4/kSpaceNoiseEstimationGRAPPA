function [gFactorMC,gFactorK,gFactorX] = ScriptRealPhantom_VD(SNR,PLOT,filepath)
%ScriptRealPhantom_VD Function to generate the figure in the paper
%comparing the g-factors for the Variable Density case with an image of
%size 128x128
% Input parameters
%   - SNR: the phantom was acquired with SNR={5,8}
%   - PLOT: flag to plot the results within the function
%   - filepath: where to save the results
%
% Output paremeters:
%   - gFactorMC: g-factor maps computed from Monte Carlo simulations
%   - gFactorK: g-factor maps computed from the k-space approach
%   - gFactorX: g-factor maps computed from the image-space approach
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

%We get the input parameters
if (nargin<1)
    SNR=5;
end
if(nargin<2)
    PLOT=0;
end

%We import the acquired water phantom
load(['../Data/Real/PhantomWaterBall_SNR' num2str(SNR) '_B0DriftRemoved.mat'])

%We define the parameters we need
Size=size(K,1);
L=size(K,3);
NEX=size(K,4);
method=1;

%We define the GRAPPA parameters
kernel=[2,3];
Z = zeros(Size,1);
Z(2:4:14) = 1;
Z(116:4:128)=1; %R=4
Z(14:3:32) = 1;
Z(98:3:116)=1; %R=3
Z(32:2:48) = 1;
Z(82:2:98)=1; %R=2
Z(48:82)=1; %R=1
LinesACS=48:82;
Lines = find(Z);
mssLines = find(1-Z);
kernelstr=[num2str(kernel(1)) 'x' num2str(kernel(2))];
N=[];
Reff=Size/length(Lines);

%We estimate the noise matrix from k-space, only considering the points
%without signal
maxValue=max(reshape(abs(K(:,:,1,1)),[],1));
find(abs(K(:,:,1,1))<0.005*maxValue);
Ind=find(abs(K(:,:,1,1))<0.005*maxValue);
KCell=num2cell(permute(K,[4 3 1 2]),[1 2]);
GammaKCell=reshape(cellfun(@(x) GammaCov(x),KCell,'UniformOutput',0),[1 1 Size*Size]);
GammaK=cell2mat(GammaKCell);
GammaK=sum(GammaK(:,:,Ind),3)/(Size*Size);
CKCell=reshape(cellfun(@(x) CCov(x),KCell,'UniformOutput',0),[1 1 Size*Size]);
CK=cell2mat(CKCell);
CK=sum(CK(:,:,Ind),3)/(Size*Size);
clear GammaKCell CKCell

%We subsampled the images and reconstruct them using GRAPPA as an
%interpolation
k_SN=zeros(Size,Size,L);
for Exi=1:NEX   
    k_SN(Lines,:,:)=K(Lines,:,:,Exi); 
    [Fc,~,N,LnsKnl,mssKnl]=recongrappa_multik([Size,Size,L],k_SN(Lines,:,:),Lines,'kernel',kernelstr,'N',N);
    Sn(:,:,:,Exi)=k2x(Fc,method);
end

%Now we estimate the correlation of noise and signal for every pixel. 
Sn_signal=mean(Sn,4);
Sn_signalCell=num2cell(permute(Sn_signal,[3 1 2]),[1]);
RsCell=cellfun(@(x) x*x',Sn_signalCell,'UniformOutput',0);

Sn_noise=Sn-repmat(Sn_signal,[1 1 1 NEX]);
Sn_noiseCell=num2cell(permute(Sn_noise,[4 3 1 2]),[1 2]);
RnCell=permute(cellfun(@(x) GammaCov(x),Sn_noiseCell,'UniformOutput',0),[1 3 4 2]);

%And with that we get the coil combination vector proposed by Walsh
WalshCell=cellfun(@(x,y) WalshFilter(x,y),RsCell,RnCell,'UniformOutput',0);
Walsh=permute(cell2mat(WalshCell),[2 3 1]);

%Now we combine the coils and extract the noise maps
ICoilsfull=k2x(K,method);
for Exi=1:NEX   
    %Coil combination
    [Ir(:,:,Exi)]=coilCombine_SingleSlice(Sn(:,:,:,Exi),Walsh);
    [Ifull(:,:,Exi)]=coilCombine_SingleSlice(ICoilsfull(:,:,:,Exi),Walsh);
end

%We estimate the noise in the fully sampled image
WalshAux=permute(Walsh,[3 1 2]);
GammaXFullCoils=GammaK/(Size^2);
CXFullCoils = CXfromK( CK,Size,L,1:Size );
CXCoilsCell=permute(num2cell(permute(CXFullCoils,[3 4 1 2]),[1 2]),[3 4 1 2]);
GammaFull=cellfun(@(y) y'*GammaXFullCoils*y,permute(num2cell(WalshAux,[1]),[2 3 1]));
CFull=cellfun(@(x,y) y'*x*conj(y),CXCoilsCell,permute(num2cell(WalshAux,[1]),[2 3 1]));
NoiseMapFull=complex(real(GammaFull+CFull),real(GammaFull-CFull))/2;
Std_Est_Full_Single=complex(sqrt(real(NoiseMapFull)),sqrt(imag(NoiseMapFull)));

%We compute the image statistics
Std_MCAcc=complex(std(real(Ir),0,3),std(imag(Ir),0,3));
Std_MCFull=complex(std(real(Ifull),0,3),std(imag(Ifull),0,3));
gFactorMC=(1/sqrt(Reff))*(Std_MCAcc./Std_MCFull);
IrCell=permute(num2cell(permute(Ir,[3 1 2]),1),[2 3 1]);
NoiseMapCross_MC=cellfun(@(x) myCrossCov(x),IrCell);
Std_MC_Coils=complex(std(real(Sn),0,4),std(imag(Sn),0,4));

%K-space estimation
[NoiseMapEstK]=NoiseGRAPPA_KSpace(GammaK,CK,kernel,L,N,Walsh,Lines,mssLines,LnsKnl,mssKnl,Size,method);
Std_Est_KSingle=complex(sqrt(real(NoiseMapEstK)),sqrt(imag(NoiseMapEstK)));
gFactorK=(1/sqrt(Reff))*(Std_Est_KSingle./Std_Est_Full_Single);

%X-space estimation
[NoiseMapEstX]=NoiseGRAPPA_XSpace(GammaK,CK,kernel,L,N,Walsh,Size,Lines,1,LinesACS);
Std_Est_XSingle=complex(sqrt(real(NoiseMapEstX)),sqrt(imag(NoiseMapEstX)));
gFactorX=(1/sqrt(Reff))*(Std_Est_XSingle./Std_Est_Full_Single);
    
%We plot the MonteCarlo and estimated maps
if (PLOT)
    range=[min(abs([gFactorMC(:);gFactorK(:);gFactorX(:)])) max(abs([gFactorMC(:);gFactorK(:);gFactorX(:)]))];
    figure
    subplot(1,3,1)
    imshow(abs(gFactorMC),[range]),colormap 'jet',colorbar,title('g-factor Map Monte Carlo');
    subplot(1,3,2)
    imshow(abs(gFactorK),[range]),colormap 'jet',colorbar,title('g-factor Map Estimated k-space');
    subplot(1,3,3)
    imshow(abs(gFactorX),[range]),colormap 'jet',colorbar,title('g-factor Map Estimated k-space');
end

%We save the results
save([filepath '/RealPhantom/Figure_VD/ResultsVD_SNR' num2str(SNR) '.mat'],'GammaK','CK','kernel','N','Walsh','Lines','mssLines','LnsKnl','mssKnl','Size','method','Reff','LinesACS','Std_MCAcc','Std_Est_KSingle','Std_Est_XSingle','Std_Est_Full_Single','Std_MCFull','gFactorMC','gFactorK','gFactorX')
end

