% Script to plot all the figures in the results section of the paper 
% [1] Rabanillo I. Aja-Fernández S. Alberola-López C., Hernando D. (2016), 
% “Exact Calculation of Noise Maps and g-Factor in GRAPPA using a k–space 
% Analysis" 

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

%We add the code to the path
addpath(genpath('../../kSpaceNoiseEstimationGRAPPA'))

%First we get the results for the synthetic phantom
AVAILABLE=1;

%Synthetic experiments figures
filepath='../Results';
if (AVAILABLE)
    load([filepath '/Synthetic/Figure_noACS/ResultsnoACS_R2_kernel2x3_SNR30_Size128.mat'])
    gFactorMC_noACSSynthetic=gFactorMC;
    gFactorK_noACSSynthetic=gFactorK;
    gFactorX_noACSSynthetic=gFactorX;

    load([filepath '/Synthetic/Figure_ACS/ResultsACS_R3_kernel2x3_SNR30_Size128.mat'])
    gFactorMC_ACSSynthetic=gFactorMC;
    gFactorK_ACSSynthetic=gFactorK;
    gFactorX_ACSSynthetic=gFactorX;    
    
    load([filepath '/Synthetic/Figure_ACS/ResultsACS_R3_kernel4x3_SNR30_Size128.mat'])
    gFactorMC_ACSSynthetic2=gFactorMC;
    gFactorK_ACSSynthetic2=gFactorK;
    gFactorX_ACSSynthetic2=gFactorX;
    
    load([filepath '/Synthetic/Figure_ACS/ResultsACS_R4_kernel2x3_SNR30_Size128.mat'])
    gFactorMC_ACSSynthetic3=gFactorMC;
    gFactorK_ACSSynthetic3=gFactorK;
    gFactorX_ACSSynthetic3=gFactorX;    
    
    load([filepath '/Synthetic/Figure_VD/ResultsVD_SNR30_Size128.mat'])
    gFactorMC_VDSynthetic=gFactorMC;
    gFactorK_VDSynthetic=gFactorK;
    gFactorX_VDSynthetic=gFactorX;    
else    
    [gFactorMC_noACSSynthetic,gFactorK_noACSSynthetic,gFactorX_noACSSynthetic]=ScriptFigure_noACS([2 3],2,30,0.1,8,128,4000,0,filepath);
    [gFactorMC_ACSSynthetic,gFactorK_ACSSynthetic,gFactorX_ACSSynthetic]=ScriptFigure_ACS([2 3],3,30,0.1,8,128,4000,0,filepath);
    [gFactorMC_ACSSynthetic2,gFactorK_ACSSynthetic2,gFactorX_ACSSynthetic2]=ScriptFigure_ACS([4 3],3,30,0.1,8,128,4000,0,filepath);
    [gFactorMC_ACSSynthetic3,gFactorK_ACSSynthetic3,gFactorX_ACSSynthetic3]=ScriptFigure_ACS([2 3],4,30,0.1,8,128,4000,0,filepath);
    [gFactorMC_VDSynthetic,gFactorK_VDSynthetic,gFactorX_VDSynthetic]=ScriptFigure_VD_Size128(30,0.1,8,4000,0,filepath);
end

figure('units','normalized','outerposition',[0 0 1 1]);
range=[0 max([abs(gFactorMC_noACSSynthetic(:));abs(gFactorK_noACSSynthetic(:));abs(gFactorX_noACSSynthetic(:))])];
subplot(5,4,1)
imshow(abs(gFactorMC_noACSSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,2)
imshow(abs(gFactorK_noACSSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,3)
imshow(abs(gFactorX_noACSSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,4)
imshow(abs(gFactorX_noACSSynthetic-gFactorK_noACSSynthetic),[0.1*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_ACSSynthetic(:));abs(gFactorK_ACSSynthetic(:));abs(gFactorX_ACSSynthetic(:))])];
subplot(5,4,5)
imshow(abs(gFactorMC_ACSSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,6)
imshow(abs(gFactorK_ACSSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,7)
imshow(abs(gFactorX_ACSSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,8)
imshow(abs(gFactorX_ACSSynthetic-gFactorK_ACSSynthetic),[0.1*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_ACSSynthetic2(:));abs(gFactorK_ACSSynthetic2(:));abs(gFactorX_ACSSynthetic2(:))])];
subplot(5,4,9)
imshow(abs(gFactorMC_ACSSynthetic2),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,10)
imshow(abs(gFactorK_ACSSynthetic2),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,11)
imshow(abs(gFactorX_ACSSynthetic2),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,12)
imshow(abs(gFactorX_ACSSynthetic2-gFactorK_ACSSynthetic2),[0.1*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_ACSSynthetic3(:));abs(gFactorK_ACSSynthetic3(:));abs(gFactorX_ACSSynthetic3(:))])];
subplot(5,4,13)
imshow(abs(gFactorMC_ACSSynthetic3),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,14)
imshow(abs(gFactorK_ACSSynthetic3),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,15)
imshow(abs(gFactorX_ACSSynthetic3),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,16)
imshow(abs(gFactorX_ACSSynthetic3-gFactorK_ACSSynthetic3),[0.1*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_VDSynthetic(:));abs(gFactorK_VDSynthetic(:));abs(gFactorX_VDSynthetic(:))])];
subplot(5,4,17)
imshow(abs(gFactorMC_VDSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,18)
imshow(abs(gFactorK_VDSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,19)
imshow(abs(gFactorX_VDSynthetic),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,20)
imshow(abs(gFactorX_VDSynthetic-gFactorK_VDSynthetic),[0.1*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');

%Real Phantom experiments figures
if (AVAILABLE)
    load([filepath '/RealPhantom/Figure_noACS/ResultsnoACS_R2_kernel2x3_SNR5.mat'])
    gFactorMC_noACSRealPhantom=gFactorMC;
    gFactorK_noACSRealPhantom=gFactorK;
    gFactorX_noACSRealPhantom=gFactorX;

    load([filepath '/RealPhantom/Figure_ACS/ResultsACS_R3_kernel2x3_SNR5.mat'])
    gFactorMC_ACSRealPhantom=gFactorMC;
    gFactorK_ACSRealPhantom=gFactorK;
    gFactorX_ACSRealPhantom=gFactorX;    
    
    load([filepath '/RealPhantom/Figure_ACS/ResultsACS_R3_kernel4x3_SNR5.mat'])
    gFactorMC_ACSRealPhantom2=gFactorMC;
    gFactorK_ACSRealPhantom2=gFactorK;
    gFactorX_ACSRealPhantom2=gFactorX;        
    
    load([filepath '/RealPhantom/Figure_ACS/ResultsACS_R4_kernel2x3_SNR5.mat'])
    gFactorMC_ACSRealPhantom3=gFactorMC;
    gFactorK_ACSRealPhantom3=gFactorK;
    gFactorX_ACSRealPhantom3=gFactorX;    
    
    load([filepath '/RealPhantom/Figure_VD/ResultsVD_SNR5.mat'])
    gFactorMC_VDRealPhantom=gFactorMC;
    gFactorK_VDRealPhantom=gFactorK;
    gFactorX_VDRealPhantom=gFactorX;    
else    
    [gFactorMC_noACSRealPhantom,gFactorK_noACSRealPhantom,gFactorX_noACSRealPhantom]=ScriptRealPhantom_noACS([2 3],2,5,0,filepath);
    [gFactorMC_ACSRealPhantom,gFactorK_ACSRealPhantom,gFactorX_ACSRealPhantom]=ScriptRealPhantom_ACS([2 3],3,5,0,filepath);
    [gFactorMC_ACSRealPhantom2,gFactorK_ACSRealPhantom2,gFactorX_ACSRealPhantom2]=ScriptRealPhantom_ACS([4 3],3,5,0,filepath);
    [gFactorMC_ACSRealPhantom3,gFactorK_ACSRealPhantom3,gFactorX_ACSRealPhantom3]=ScriptRealPhantom_ACS([2 3],4,5,0,filepath);
    [gFactorMC_VDRealPhantom,gFactorK_VDRealPhantom,gFactorX_VDRealPhantom]=ScriptRealPhantom_VD(5,0,filepath);
end

figure('units','normalized','outerposition',[0 0 1 1]);
range=[0 max([abs(gFactorMC_noACSRealPhantom(:));abs(gFactorK_noACSRealPhantom(:));abs(gFactorX_noACSRealPhantom(:))])];
subplot(5,4,1)
imshow(abs(gFactorMC_noACSRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,2)
imshow(abs(gFactorK_noACSRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,3)
imshow(abs(gFactorX_noACSRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,4)
imshow(abs(gFactorX_noACSRealPhantom-gFactorK_noACSRealPhantom),[0.05*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_ACSRealPhantom(:));abs(gFactorK_ACSRealPhantom(:));abs(gFactorX_ACSRealPhantom(:))])];
subplot(5,4,5)
imshow(abs(gFactorMC_ACSRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,6)
imshow(abs(gFactorK_ACSRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,7)
imshow(abs(gFactorX_ACSRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,8)
imshow(abs(gFactorX_ACSRealPhantom-gFactorK_ACSRealPhantom),[0.05*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_ACSRealPhantom2(:));abs(gFactorK_ACSRealPhantom2(:));abs(gFactorX_ACSRealPhantom2(:))])];
subplot(5,4,9)
imshow(abs(gFactorMC_ACSRealPhantom2),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,10)
imshow(abs(gFactorK_ACSRealPhantom2),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,11)
imshow(abs(gFactorX_ACSRealPhantom2),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,12)
imshow(abs(gFactorX_ACSRealPhantom2-gFactorK_ACSRealPhantom2),[0.05*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_ACSRealPhantom3(:));abs(gFactorK_ACSRealPhantom3(:));abs(gFactorX_ACSRealPhantom3(:))])];
subplot(5,4,13)
imshow(abs(gFactorMC_ACSRealPhantom3),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,14)
imshow(abs(gFactorK_ACSRealPhantom3),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,15)
imshow(abs(gFactorX_ACSRealPhantom3),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,16)
imshow(abs(gFactorX_ACSRealPhantom3-gFactorK_ACSRealPhantom3),[0.05*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');
range=[0 max([abs(gFactorMC_VDRealPhantom(:));abs(gFactorK_VDRealPhantom(:));abs(gFactorX_VDRealPhantom(:))])];
subplot(5,4,17)
imshow(abs(gFactorMC_VDRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Monte Carlo');
subplot(5,4,18)
imshow(abs(gFactorK_VDRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Estimated k-space');
subplot(5,4,19)
imshow(abs(gFactorX_VDRealPhantom),[range]),colormap 'jet',colorbar,title('g-factor Estimated Image-space');
subplot(5,4,20)
imshow(abs(gFactorX_VDRealPhantom-gFactorK_VDRealPhantom),[0.05*range]),colormap 'jet',colorbar,title('g-factor Differences k-space and Image-space');