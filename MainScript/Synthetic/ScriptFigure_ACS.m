function [gFactorMC,gFactorK,gFactorX]=ScriptFigure_ACS(kernel,R,SNR,rho,L,Size,NEX,PLOT,filepath)
%ScriptFigure_noACS Function to generate the figure in the paper
%comparing the g-factors in the case when the ACS region is present
%   - kernel: [kx,ky] vector containing the size of the kernel
%   - R: acceleration factor
%   - SNR: signal to noise ratio in the original fully-sampled images
%   - rho: correlation coefficient accross channels
%   - L: number of coils
%   - Size: size of the fully-sampled image
%   - NEX: number or repetitions for the Monte Carlo simulations
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
if(nargin<1)
    kernel=[2 3];
end
if (nargin<2)
    R=2;
end
if (nargin<3)
    SNR=25;
end
if (nargin<4)
    rho=0.1;
end
if (nargin<5)
    L=8;
end
if (nargin<6)
    Size=64;
end
if (nargin<7)
    NEX=4000;
end
if (nargin<8)
    PLOT=0;
end

%This is to compare both estimations against the Monte-Carlo maps. We
%define the sensitivity maps and obtain the fully-sampled coil images
% load('/home/inaki/Escritorio/ProyectoMadison/Data_LiverIron/NoiseGRAPPA/TestBrainWeb/RealBrain256.mat');
load('../Data/Synthetic/BrainWeb256.mat');
I=I(256/Size:256/Size:end,256/Size:256/Size:end);
MapW=sensitivity_map([Size,Size],L); %Sensitivity Map
It=repmat(double(I),[1,1,L]).*MapW;

%We define the stationary noise matrices in the k-space
sigmaK=(mean(abs(I(:)))/SNR)*Size;
VxxK=sigmaK^2*(eye(L)+(ones(L)-eye(L))*rho);
VyyK=sigmaK^2*(eye(L)+(ones(L)-eye(L))*rho);
MAux=ones(L);
IndTriLowDiag=find(tril(MAux))';
IndTriLow=find(tril(MAux,-1))';
IndTriUp=setdiff(1:L^2,IndTriLowDiag); 
MAux2=zeros(L);
MAux2(IndTriUp)=IndTriUp;
MAux2=MAux2';
MAux2=zeros(L);
MAux2(IndTriUp)=IndTriUp;
MAux2=MAux2';
IndTriUp=MAux2(find(MAux2))';
VxyK=zeros(L);
VxyK(IndTriLow)=2*(randi(2,length(IndTriLow),1)-1.5)*rho;
VxyK(IndTriUp)=-VxyK(IndTriLow);
VxyK=sigmaK^2*VxyK;
VyxK=VxyK.';
VTot=[VxxK VxyK;VyxK VyyK];
GammaK=complex(VxxK+VyyK,VyxK-VxyK);
CK=complex(VxxK-VyyK,VyxK+VxyK);
W=chol(VTot)';

%We define the GRAPPA parameters
Z = zeros(Size,1);
LinesSub= 2:R:Size;
Z( LinesSub) = 1;  
Nacs=32;
LinesACS=(Size/2) + (-(Nacs/2-2):Nacs/2);
LinesACSSub=LinesACS(ismember(LinesACS,LinesSub));
LinesACS=LinesACSSub(1):LinesACSSub(end);
Z( LinesACS ) = 1;	%ACS lines 
Lines = find(Z);
Reff=Size/length(Lines);
mssLines = find(1-Z);
kernelstr=[num2str(kernel(1)) 'x' num2str(kernel(2))];
method=1;
N=[];
Walsh=[];

%We define the matrices to compute the Monte Carlo results
NoiseMeanAcc=zeros([Size Size]);
NoiseMSEAcc=zeros([Size Size]);
NoiseMSECrossAcc=zeros([Size Size]);
NoiseMeanFull=zeros([Size Size]);
NoiseMSEFull=zeros([Size Size]);

%We transform the image to k-space
K=x2k(It);
k_SN=zeros(Size,Size,L);
for Exi=1:NEX
    %We corrupt the data with synthetic noise
    Nc=randn([2*L,Size^2]);
    Nc2=W*Nc;
    Nc3=complex(Nc2(1:L,:),Nc2((L+1):(2*L),:));
    Nc4=permute(reshape(Nc3,[L,Size,Size]),[2 3 1]);
    kN=K+Nc4;
    
    %We reconstruct the image with no acceleration
    Snfull=k2x(kN,method);
    [Ifull,Walsh]=coilCombine_SingleSlice(Snfull,Walsh);
    NoiseMeanFull=NoiseMeanFull+Ifull;
    NoiseMSEFull=NoiseMSEFull+complex(real(Ifull).^2,imag(Ifull).^2);
    
    %We reconstruct the subsampled image with GRAPPA as a convolution
    k_SN(Lines,:,:)=kN(Lines,:,:); 
    [Fc,~,N,LnsKnl,mssKnl]=recongrappa_multik([Size,Size,L],k_SN(Lines,:,:),Lines,'kernel',kernelstr,'dks,',R,'N',N);
    Sn=k2x(Fc,method);
    [Ir,Walsh]=coilCombine_SingleSlice(Sn,Walsh); %Deterministic
    NoiseMeanAcc=NoiseMeanAcc+Ir;
    NoiseMSEAcc=NoiseMSEAcc+complex(real(Ir).^2,imag(Ir).^2);
    NoiseMSECrossAcc=NoiseMSECrossAcc+real(Ir).*imag(Ir);    
end
clear Nc Nc2 k_SN Fc Sn Ir

%We get the Monte Carlo noise maps
NoiseMeanAcc=NoiseMeanAcc/NEX;
NoiseMSEAcc=NoiseMSEAcc/NEX;
NoiseMSECrossAcc=NoiseMSECrossAcc/NEX;
NoiseMapMCAcc=NoiseMSEAcc-complex(real(NoiseMeanAcc).^2,imag(NoiseMeanAcc).^2);
NoiseMapCrossAcc=NoiseMSECrossAcc-real(NoiseMeanAcc).*imag(NoiseMeanAcc);
Std_MCAcc=complex(sqrt(real(NoiseMapMCAcc)),sqrt(imag(NoiseMapMCAcc)));
NoiseMeanFull=NoiseMeanFull/NEX;
NoiseMSEFull=NoiseMSEFull/NEX;
NoiseMapMCFull=NoiseMSEFull-complex(real(NoiseMeanFull).^2,imag(NoiseMeanFull).^2);
Std_MCFull=complex(sqrt(real(NoiseMapMCFull)),sqrt(imag(NoiseMapMCFull)));
gFactorMC=(1/sqrt(Reff))*(Std_MCAcc./Std_MCFull);

%We estimate the noise in the fully sampled image
WalshAux=permute(Walsh,[3 1 2]);
GammaXFullCoils=GammaK/(Size^2);
CXFullCoils = CXfromK( CK,Size,L,1:Size );
CXCoilsCell=permute(num2cell(permute(CXFullCoils,[3 4 1 2]),[1 2]),[3 4 1 2]);
GammaFull=cellfun(@(y) y'*GammaXFullCoils*y,permute(num2cell(WalshAux,[1]),[2 3 1]));
CFull=cellfun(@(x,y) y'*x*conj(y),CXCoilsCell,permute(num2cell(WalshAux,[1]),[2 3 1]));
NoiseMapFull=complex(real(GammaFull+CFull),real(GammaFull-CFull))/2;
Std_Est_Full=complex(sqrt(real(NoiseMapFull)),sqrt(imag(NoiseMapFull)));

%K-space estimation
[NoiseMapK,NoiseMapCrossK]=NoiseGRAPPA_KSpace(GammaK,CK,kernel,L,N,Walsh,Lines,mssLines,LnsKnl,mssKnl,Size,method);
Std_Est_K=complex(sqrt(real(NoiseMapK)),sqrt(imag(NoiseMapK)));
gFactorK=(1/sqrt(Reff))*(Std_Est_K./Std_Est_Full);

%X-space estimation
[NoiseMapX,NoiseMapCrossX]=NoiseGRAPPA_XSpace(GammaK,CK,kernel,L,N,Walsh,Size,Lines,1,LinesACS);
Std_Est_X=complex(sqrt(real(NoiseMapX)),sqrt(imag(NoiseMapX)));
gFactorX=(1/sqrt(Reff))*(Std_Est_X./Std_Est_Full);

%We plot the MonteCarlo and estimated maps
if (PLOT)
    range=[min([abs(NoiseMapMC(:));abs(NoiseMapX(:));abs(NoiseMapK(:))]) max([abs(NoiseMapMC(:));abs(NoiseMapX(:));abs(NoiseMapK(:))])];
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1)
    imshow(abs(NoiseMapMC),[range]),colormap 'jet',colorbar,title('Noise Map Monte Carlo');
    subplot(1,3,2)
    imshow(abs(NoiseMapK),[range]),colormap 'jet',colorbar,title('Noise Map Estimated k-space');
    subplot(1,3,3)
    imshow(abs(NoiseMapX),[range]),colormap 'jet',colorbar,title('Noise Map Estimated x-space');
    print(['/home/inaki/Escritorio/ProyectoMadison/CircularConv/Figure_ACS/NoiseMap_ACS_R' num2str(R) '_kernel' kernelstr '_SNR' num2str(SNR) '_Size' num2str(Size)],'-dpng')

    range=[min([abs(gFactorMC(:));abs(gFactorX(:));abs(gFactorK(:))]) max([abs(gFactorMC(:));abs(gFactorX(:));abs(gFactorK(:))])];
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1)
    imshow(abs(gFactorMC),[range]),colormap 'jet',colorbar,title('g-factor Map Monte Carlo');
    subplot(1,3,2)
    imshow(abs(gFactorK),[range]),colormap 'jet',colorbar,title('g-factor Map Estimated k-space');
    subplot(1,3,3)
    imshow(abs(gFactorX),[range]),colormap 'jet',colorbar,title('g-factor Map Estimated x-space');
    print(['/home/inaki/Escritorio/ProyectoMadison/CircularConv/Figure_ACS/g-factor_ACS_R' num2str(R) '_kernel' kernelstr '_SNR' num2str(SNR) '_Size' num2str(Size)],'-dpng')

    range=[min([NoiseMapCrossAcc(:);NoiseMapCrossK(:);NoiseMapCrossX(:)]) max([NoiseMapCrossAcc(:);NoiseMapCrossK(:);NoiseMapCrossX(:)])]
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1)
    imshow(abs(NoiseMapCrossAcc),[range]),colormap 'jet',colorbar,title('Noise Map Cross Monte Carlo');
    subplot(1,3,2)
    imshow(abs(NoiseMapCrossK),[range]),colormap 'jet',colorbar,title('Noise Map Cross Estimated k-space');
    subplot(1,3,3)
    imshow(abs(NoiseMapCrossX),[range]),colormap 'jet',colorbar,title('Noise Map Cross Estimated x-space');
    print([filepath '/Figure_noACS/NoiseMapCross_ACS_ACS_R' num2str(R) '_kernel' kernelstr '_SNR' num2str(SNR) '_Size' num2str(Size)],'-dpng')
end

%We save the results
save([filepath '/Synthetic/Figure_ACS/ResultsACS_R' num2str(R) '_kernel' kernelstr '_SNR' num2str(SNR) '_Size' num2str(Size) '.mat'],'GammaK','CK','kernel','SNR','L','N','Walsh','Lines','mssLines','LnsKnl','mssKnl','Size','method','R','Reff','LinesACS','NoiseMapMCAcc','NoiseMapK','NoiseMapX','NoiseMapFull','NoiseMapMCFull','gFactorMC','gFactorK','gFactorX','NoiseMapCrossAcc','NoiseMapCrossK','NoiseMapCrossX')
end