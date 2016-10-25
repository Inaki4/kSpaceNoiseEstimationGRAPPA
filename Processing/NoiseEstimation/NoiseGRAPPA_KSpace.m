function [NoiseMap,NoiseMapXY,NoiseMapCoils] = NoiseGRAPPA_KSpace(Gamma,C,kernel,L,N,Walsh,acqLines,mssLines,LnsKnl,mssKnl,Size,method)
% function [NoiseMap,ColCorrGamma,ColCorrGammaFFT,RowCorr,RowCorrGammaFFT,CoilCrrGamma] = NoiseGRAPPA(Gamma,kernel,R,L,N,Walsh,acqLines,mssLines,LnsKnl,mssKnl,Size,method)
%NOISEGRAPPA Generates a noise map based on the GRAPPA reconstruction
%performed in the k-space. 
%Input parameters:
%   - Gamma,C: Covariance matrices for a complex multivariate gaussian (see:
%   https://en.wikipedia.org/wiki/Complex_normal_distribution). 
%   - kernel=[ky kx]: size of the kernel used for the reconstruction.
%   - L: number of coils
%   - N: kernels used for the reconstruction.
%   - Walsh: 3D array [Size,Size,L] with the vector to combine all the
%   coils for each pixel in order to get the final image.
%   - acqLines: vector with the indices of the acquired lines from which 
%   the reconstruction is performed.
%   - mssLines: vector with the indices of the missing lines to fill when 
%   the reconstruction is performed
%   - LnsKnl: cell array of length(mssLines) containing for each entry the 
%   indices of the lines from which the current missing line is
%   reconstructed.
%   - mssKnl: cell array of length(mssLines) containing for each entry the 
%   index of the kernel used to reconstruct the current missing line.
%   - Size of the image.
%   - method: way to compute the IDFFT. Default is method=1
%
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

nkx=kernel(2);
nky=kernel(1);
Qk=nkx*nky;
M=length(N);
D=(2*nkx-1)*(2*nky-1);

%==========================================================================
%------------------------BUILD THE KERNEL MATRICES ------------------------
%==========================================================================
%First we build the kernel matrices. They are LxL matrices where there
%i-row contains the weights to build the i-coil from the values for the
%analyzed pixel along the coil dimension.
%We fill the non-shifted kernel
W=zeros(L,L,Qk,D*M); %The order in the forth dimension is we first place all the shifted versions of each kernel
WGamma=zeros(L,L,Qk,M);
WC=zeros(L,L,Qk,M);
IndQ=cell(D,1); %Stores the linear index of the pixels in the kernel preserved when shifting.
for m=1:M
    Waux=zeros(nkx*L,L,nky);
    for y=1:nky
        Waux(:,:,y)=N(m).n((1:nkx*L)+(y-1)*nkx*L,:);
    end
    Waux=permute(Waux,[3 1 2]);
    Waux=reshape(Waux,[nky,nkx,L,L]);
    
    %Now we fill all the shifted kernels
    for d=[1:D]
        Waux2=zeros(nky,nkx,L,L);
        [dy,dx]=ind2sub([2*nky-1,2*nkx-1],d);
        dm=sub2ind([D,M],d,m);
        my=max(dy-(nky-1),1);
        My=min(dy,nky);
        mx=max(dx-(nkx-1),1);
        Mx=min(dx,nkx);
        
        if m==1
            [AuxY,AuxX]=meshgrid(my:My,mx:Mx);
            IndQ{dm}=sub2ind([nky,nkx],AuxY(:),AuxX(:));
            clear AuxY AuxX;
        end
        
        Waux2(my:My,mx:Mx,:,:)=Waux((nky+1-My):(nky+1-my),(nkx+1-Mx):(nkx+1-mx),:,:);
        Waux2=reshape(Waux2,[Qk,L,L]);
        Waux2=permute(Waux2,[2 3 1]); %[L,L,Qk] array with the weights for all coils to reconstruct each coil for each pixel in the kernel m with a shift [dx,dy]
        W(:,:,:,dm)=Waux2;
    end
    clear my My mx Mx Waux2 
    Waux=reshape(Waux,[Qk,L,L]);
    Waux=permute(Waux,[2 3 1]);
    WGammaaux=cellfun(@(x) (x.')*Gamma,num2cell(Waux,[1 2]),'UniformOutput',false);
    WCaux=cellfun(@(x) (x.')*C,num2cell(Waux,[1 2]),'UniformOutput',false);
    WGamma(:,:,:,m)=cat(3,WGammaaux{:});
    WC(:,:,:,m)=cat(3,WCaux{:});
    clear Waux WGammaaux
end
clear d m dm dy dx mx Mx my My

%==========================================================================
%------------------------BUILD COLUMN CORRELATIONS-------------------------
%==========================================================================
IndFilledLines=find(mssKnl);
FilledLines=reshape(mssLines(IndFilledLines),1,[]);
MattLnsKnl=cell2mat(LnsKnl(IndFilledLines));
cnt1=0;
NLnsCorr=zeros(Size,1);

%The reference point is always the unshifted one, whose index is
WGammaW_Data=cell(M,M*D);
WCW_Data=cell(M,M*D);

%We build the index of the lower half of the correlation matrix that we are
%going to compute due to its symmetry
[IndX,IndY]=meshgrid(1:L,1:L);
IndCorr=find(IndX<=IndY);
clear IndX IndY
PntCrr=length(IndCorr);
ColCorrGamma=zeros(Size,Size*(2*nkx-1),PntCrr);
ColCorrC=zeros(Size,Size*(2*nkx-1),PntCrr);

%We define the cell that will contain the kernels of the correlated lines
%with each of the reference line
KnlsRefLn=cell(1,length(FilledLines));
IndxRefLnCrrLns=cell(1,length(FilledLines));

%We define the lines used to reconstruct
IndACS=find(diff(acqLines)==1);
RecLines=acqLines;
if ~isempty(IndACS)
    AddedLnKnl=[];
    if ((kernel(1)/2-1)>0)
        Rinf=acqLines(IndACS(1))-acqLines(IndACS(1)-1);
        Rsup=acqLines(IndACS(end)+2)-acqLines(IndACS(end)+1);
        AddedLnKnl=[acqLines(IndACS(1))+Rinf*(1:(kernel(1)/2-1))';acqLines(IndACS(end)+1)-Rsup*((kernel(1)/2-1):-1:1)'];
    end
    RecLines=[acqLines(1:(IndACS(1))); AddedLnKnl; acqLines((1+IndACS(end)):end)];
end
clear IndACS AddedLnKnl Rinf Rsup AddedLnKnl

%Now we loop thourg the reconstructed lines
for refLn=FilledLines
    cnt1=cnt1+1;
    
    %We find the lines the current line correlates with
    mssLnsCorr=FilledLines(find(sum(ismember(MattLnsKnl,MattLnsKnl(cnt1,:)),2)>0)); %mssLines((fFilledLine-1)+find(sum(ismember(MattLnsKnl,MattLnsKnl(cnt1,:)),2)>0));
    acqLnsCorr=MattLnsKnl(cnt1,:);
    
    %We get the current line kernel
    m=mssKnl(find(mssLines==refLn)); %=mssKnl((fFilledLine-1)+cnt1); 
    
    %We get the kernel of the missing lines the current one correlates with
    [~,IndmssLnsCorr]=intersect(mssLines,mssLnsCorr);
    IndxLn=find(mssLnsCorr==refLn);
    KnlmssLnsCorr=mssKnl(IndmssLnsCorr);
    MatchLns=cellfun(@(x,y) isequal(x,[m; KnlmssLnsCorr])&&isequal(y,IndxLn),KnlsRefLn,IndxRefLnCrrLns);
    MatchLn=FilledLines(find(MatchLns,1));
    KnlsRefLn{cnt1}=[m;KnlmssLnsCorr];
    IndxRefLnCrrLns{cnt1}=IndxLn; %This contains the position of the reference line with respect to the lines it correlates with
    
    %If the kernel of the current line and all the lines it correlates with
    %matches those of a different reference line, the correlation vector is
    %just a shifted version of it, so we don't have to compute it AGAIN. If
    %it is not the case, we need to do so. The situation has to be the
    %same, which means that the position of the reference line within the
    %lines it correlates with has to be the same.
    if isempty(MatchLn)
        for jKnlY=1:length(mssLnsCorr)
            shftLn=mssLnsCorr(jKnlY);
            acqLnsshftLn=MattLnsKnl(find(FilledLines==shftLn),:); %Acquired lines for the shifted line
            
            %We determine if the first line in the kernel for both shifted 
            %and reference lines ir reached by cycling 
            LnsInt=intersect(acqLnsCorr,acqLnsshftLn);
            LnInt=LnsInt(1); %We get one line in the kernel of both
            dRef=[acqLnsCorr(1)-LnInt LnInt-acqLnsCorr(1)];
            dShft=[acqLnsshftLn(1)-LnInt LnInt-acqLnsshftLn(1)]; %We get the distance going down and up to reach that shared line for both. A negative distance means cyclic
            [~,IndmindRefLnInt]=min(mod(dRef,Size));
            [~,IndmindShftLnInt]=min(mod(dShft,Size)); %We get the index of the minimum distance to reach the shared line
            
            %The lines are connected by cycling if the shortest way to
            %reach both of is negative in one and positive in the other. So
            %if the sign is the same, there is no cycling.
            if(mysign(dRef(IndmindRefLnInt))==mysign(dShft(IndmindShftLnInt)))
                d1stLns=find(RecLines==acqLnsshftLn(1))-find(RecLines==acqLnsCorr(1));
            else
                min1stLn=min(find(RecLines==acqLnsCorr(1)),find(RecLines==acqLnsshftLn(1)));
                RecReord=circshift(RecLines,-min1stLn);
                d1stLns=find(RecReord==acqLnsshftLn(1))-find(RecReord==acqLnsCorr(1));
            end
            sy=nky+d1stLns;
%             sy=nky+floor((acqLnsshftLn(1)-acqLnsCorr(1))/mssR((fFilledLine-1)+cnt1));%sy=nky+find(acqLines==acqLnsshftLn(1))-find(acqLines==acqLnsCorr(1)); %sy=nky+floor((acqLnsshftLn(1)-acqLnsCorr(1))/R);
            n=mssKnl(find(mssLines==shftLn));
            
            %We get the distance 
            for sx=1:(2*nkx-1)
                s=sub2ind([2*nky-1,2*nkx-1],sy,sx);
                ops=sub2ind([2*nky-1,2*nkx-1],2*nky-sy,2*nkx-sx);
                sn=sub2ind([D,M],s,n);
                opsm=sub2ind([D,M],ops,m);
                if (isempty(WGammaW_Data{m,sn}) && isempty(WGammaW_Data{n,opsm}))
                    %The correlation matrix between pixels is
                    WGammaAux=WGamma(:,:,IndQ{s},m);
                    WCAux=WC(:,:,IndQ{s},m);                    
                    WAux=W(:,:,IndQ{s},sn);
                    WGammaW=cellfun(@(x,y) x*conj(y),num2cell(WGammaAux,[1 2]),num2cell(WAux,[1 2]),'UniformOutput',false);
                    WCW=cellfun(@(x,y) x*y,num2cell(WCAux,[1 2]),num2cell(WAux,[1 2]),'UniformOutput',false);
                    WGammaW=cell2mat(WGammaW);
                    WCW=cell2mat(WCW);
                    WGammaW=sum(WGammaW,3);
                    WCW=sum(WCW,3);                    
                    WGammaW_Data{m,sn}=WGammaW;
                    WCW_Data{m,sn}=WCW;                  
                elseif (~isempty(WGammaW_Data{m,sn}))
                    WGammaW=WGammaW_Data{m,sn};
                    WCW=WCW_Data{m,sn};
                elseif (~isempty(WGammaW_Data{n,opsm}))
                    WGammaW=WGammaW_Data{n,opsm}';
                    WCW=WCW_Data{n,opsm}.';
                end
                %We fill the data
                ColCorrGamma(refLn,Size*(sx-1)+shftLn,:)=reshape(WGammaW(IndCorr),[1 1 PntCrr]);
                ColCorrC(refLn,Size*(sx-1)+shftLn,:)=reshape(WCW(IndCorr),[1 1 PntCrr]);
            end
        end
        clear jKnlY shftLn sy n sx s ops sn opsm WGammaAux WCAux WAux WGammaW WCW acqLnsshftLn
        
        %Now we fill the correlation the acquired line, which is simply the
        %kernel. We save it in the order ky,kx.
        for sy=1:nky %sy=1:length(acqLnsCorr)
            shftLn=acqLnsCorr(sy);
            for sx=1:nkx
                sysx=sub2ind([nky,nkx],sy,sx);
                WGammaaux=WGamma(:,:,sysx,m);
                WCaux=WC(:,:,sysx,m);
                ColCorrGamma(refLn,Size*((nkx-1)/2+sx-1)+shftLn,:)=reshape(WGammaaux(IndCorr),[1 1 PntCrr]);
                ColCorrC(refLn,Size*((nkx-1)/2+sx-1)+shftLn,:)=reshape(WCaux(IndCorr),[1 1 PntCrr]);
            end
        end
        clear sy shftLn sx sysx WGammaaux WCaux
    %If the current number of lines matches the previous case, we have a
    %shifted version
    else
        ColCorrGamma(refLn,:,:)=cat(2,zeros(1,refLn-MatchLn,PntCrr),ColCorrGamma(MatchLn,1:end-(refLn-MatchLn),:));
        ColCorrC(refLn,:,:)=cat(2,zeros(1,refLn-MatchLn,PntCrr),ColCorrC(MatchLn,1:end-(refLn-MatchLn),:));
    end
end
clear MatchLn MatchLns refLn KnlsRefLn KnlmssLnsCorr IndmssLnsCorr cnt1 mssLnsCorr acqLnsCorr m RecReord

%We turn it into a NxNx(2*nkx-1)*PntCrr 3D-array of blocks
ColCorrGamma=reshape(ColCorrGamma,[Size,Size,(2*nkx-1)*PntCrr]);
ColCorrC=reshape(ColCorrC,[Size,Size,(2*nkx-1)*PntCrr]);

%Now we fill the lines in the kernel, which stay the same. The symmetry for
%the kernel is respect the central column.
AuxGamma=Gamma(IndCorr);
AuxC=C(IndCorr);
AuxGamma=reshape(AuxGamma,[1,1,PntCrr]);
AuxC=reshape(AuxC,[1,1,PntCrr]);
IndOwn=nkx+(2*nkx-1)*(0:PntCrr-1);

IndKnl=(nkx-1)/2+(0:(PntCrr-1))*(2*nkx-1);
for refLn=acqLines'
    %First we fill the correlation with itself, which is Gamma
    ColCorrGamma(refLn,refLn,IndOwn)=AuxGamma;    
    ColCorrC(refLn,refLn,IndOwn)=AuxC;    
    
    %We get all the lines that are reconstructed from the current one and
    %loop on them
    IndxRecFrom=find(any(MattLnsKnl==refLn,2));
    RecFrom=FilledLines(IndxRecFrom);
    cnt1=1;
    for shftLn=RecFrom
        %We get the current line kernel
        m=mssKnl(find(mssLines==shftLn));
        sy=find(MattLnsKnl(IndxRecFrom(cnt1),:)==refLn);
        for sx=1:nkx
            sysx=sub2ind([nky,nkx],sy,sx);
            WGammaaux=WGamma(:,:,sysx,m)';
            WCaux=WC(:,:,sysx,m).';
            ColCorrGamma(refLn,shftLn,IndKnl+(nkx+1-sx))=reshape(WGammaaux(IndCorr),[1 1 PntCrr]);
            ColCorrC(refLn,shftLn,IndKnl+(nkx+1-sx))=reshape(WCaux(IndCorr),[1 1 PntCrr]);
        end
        cnt1=cnt1+1;
    end   
end
clear AuxGamma AuxC IndAux IndOwn IndxRecFrom cnt1 WGammaaux WCaux

%==========================================================================
%-----------------------GET COLUMN FFT CORRELATIONS------------------------
%==========================================================================
%We just have to transform the blocks
ColCorrGammaFFT=NoiseFFTGamma(ColCorrGamma);
ColCorrCFFT=NoiseFFTC(ColCorrC);
clear ColCorrGamma ColCorrC;

%==========================================================================
%-----------------BUILD ROW CORRELATIONS IN HYBRID SPACE-------------------
%==========================================================================
%%For each row we build its correlation matrix. We only care about the 
%correlations within the row for all the coils.
IndDiag=cumsum([1:(Size+1):Size^2; Size^2.*ones((2*nkx-1)*PntCrr-1,Size)]);%#diagonal indices
IndDiag=reshape(IndDiag,[(2*nkx-1) PntCrr*Size]);
RowCorrGamma=zeros(Size,PntCrr*Size);
RowCorrGamma(1+mod((nkx-1):-1:(-nkx+1),Size),:)=ColCorrGammaFFT(IndDiag);
RowCorrC=zeros(Size,PntCrr*Size);
RowCorrC(1+mod((nkx-1):-1:(-nkx+1),Size),:)=ColCorrCFFT(IndDiag);
RowCorrGammaFFT=ifft(RowCorrGamma,[],1);
RowCorrCFFT=ifft(RowCorrC,[],1);
clear IndDiag RowCorrGamma RowCorrC

%==========================================================================
%-------------------GET FINAL CORRELATIONS FOR EACH COIL-------------------
%==========================================================================
CoilTriangCrrGamma=permute(reshape(RowCorrGammaFFT,[Size,PntCrr,Size]),[3 1 2]);
CoilTriangCrrC=zeros(Size,Size,PntCrr);
CoilTriangCrrC(:,[1 Size/2+1],:)=permute(reshape(RowCorrCFFT([1 Size/2+1],:),[2,PntCrr,Size]),[3 1 2]);
clear RowCorrGammaFFT RowCorrCFFT

%==========================================================================

%==========================================================================
%-------------------------GET FINAL CORRELATIONS---------------------------
%==========================================================================
%The last step is to combine all coils using the Walsh method
MAux=ones(L);
IndTriLowDiag=find(tril(MAux))';
IndTriLow=find(tril(MAux,-1))';
IndTriUp=setdiff(1:L^2,IndTriLowDiag);
MAux2=zeros(L);
MAux2(IndTriUp)=IndTriUp;
MAux2=MAux2';
IndTriUp=MAux2(find(MAux2))';
clear MAux MAux2
IndTriLowDiag=repmat(IndTriLowDiag,[1 Size^2]);
IndTriLow=repmat(IndTriLow,[1 Size^2]);
IndTriUp=repmat(IndTriUp,[1 Size^2]);
IndAuxLowDiag=reshape(repmat((0:(Size^2-1))*L^2,[PntCrr 1]),1,[]);
IndAuxnonDiag=reshape(repmat((0:(Size^2-1))*L^2,[L^2-PntCrr 1]),1,[]);
IndTriLowDiag=IndTriLowDiag+IndAuxLowDiag;
IndTriLow=IndTriLow+IndAuxnonDiag;
IndTriUp=IndTriUp+IndAuxnonDiag;
AuxGamma=permute(CoilTriangCrrGamma,[3 1 2]);
AuxC=permute(CoilTriangCrrC,[3 1 2]);
CoilCrrGamma=zeros(L^2,Size,Size);
CoilCrrC=zeros(L^2,Size,Size);
CoilCrrGamma(IndTriLowDiag)=AuxGamma;
CoilCrrC(IndTriLowDiag)=AuxC;
CoilCrrGamma(IndTriUp)=CoilCrrGamma(IndTriLow)';
CoilCrrC(IndTriUp)=CoilCrrC(IndTriLow).';
CoilCrrGamma=reshape(CoilCrrGamma,[L,L,Size,Size]);
CoilCrrC=reshape(CoilCrrC,[L,L,Size,Size]);
WalshAux=permute(Walsh,[3 1 2]);
NoiseMapGamma=cellfun(@(x,y) y'*x*y,permute(num2cell(CoilCrrGamma,[1 2]),[3 4 1 2]),permute(num2cell(WalshAux,[1]),[2 3 1]));
NoiseMapC=cellfun(@(x,y) y'*x*conj(y),permute(num2cell(CoilCrrC,[1 2]),[3 4 1 2]),permute(num2cell(WalshAux,[1]),[2 3 1]));
NoiseMap=complex(real(NoiseMapGamma+NoiseMapC),real(NoiseMapGamma-NoiseMapC))/2;
NoiseMapXY=imag(-NoiseMapGamma+NoiseMapC)/2;

%We also compute the noise map for each coil, previous to combination
NoiseCoilsGamma=reshape(CoilCrrGamma,[L L Size*Size]);
NoiseCoilsC=reshape(CoilCrrC,[L L Size*Size]);
IndDiag=cumsum([1:(L+1):L^2; L^2.*ones(Size^2-1,L)])';
NoiseCoilsGamma=NoiseCoilsGamma(IndDiag);
NoiseCoilsC=NoiseCoilsC(IndDiag);
NoiseCoilsGamma=permute(reshape(NoiseCoilsGamma,[L Size Size]),[2 3 1]);
NoiseCoilsC=permute(reshape(NoiseCoilsC,[L Size Size]),[2 3 1]);
NoiseMapCoils=complex(real(NoiseCoilsGamma+NoiseCoilsC),real(NoiseCoilsGamma-NoiseCoilsC))/2;
clear IndTriLowDiag IndTriLow IndTriUp IndAuxLowDiag IndAuxnonDiag AuxGamma WalshAux
clear RowCorrGammaFFT CoilCrrGamma NoiseMapGamma NoiseMapC

%Depending on the method to order the elements, which is related to the
%output format data, we may have to rearrange them
switch method
    case 2
        NoiseMap=ifftshift(NoiseMap);
    case 3
        NoiseMap=fftshift(NoiseMap,2);
end


end

%Function to compute the noise propagation through an iFFT for matrix Gamma
function x=NoiseFFTGamma(k)
[Mk,Mx,Mz]=size(k);
for ii=1:Mz
    x(:,:,ii)=ifft(fftshift(ifft(fftshift(k(:,:,ii),1),[],1)',1),[],1)';
end
end

%Function to compute the noise propagation through an iFFT for matrix C
function x=NoiseFFTC(k)
[Mk,Mx,Mz]=size(k);
for ii=1:Mz
    x(:,:,ii)=ifft2(fftshift(fftshift(k(:,:,ii),1),2));
end
end

%Function to compute the circulant matrix from cyclical shifts of the given
%vector
function M=ShiftedCov(vcorr,d,nkx,Size)
%We get all the lines that only correlate with adjacent points
AdjCorr=convmtx(vcorr,Size-(2*nkx-2));

%We get the lines that correlate with points at the other edge of the line
vcorr_pad=[vcorr zeros(1,Size-(2*nkx-1))];
CircCorr=cell2mat(arrayfun(@(x) circshift(vcorr_pad.',x).',d,'UniformOutput',0));

%We combine all of them
M=[CircCorr(1:(nkx-1),:);AdjCorr;CircCorr(nkx:end,:)];
end

%Function to compute the sign of a number
function y=mysign(x)
    y=-1+2*(x>=0);
end