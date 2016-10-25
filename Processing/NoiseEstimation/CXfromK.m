function [ CX ] = CXfromK( CK,Size,L,Lines )
%CXFROMK Computes the C matrix (see:
%https://en.wikipedia.org/wiki/Complex_normal_distribution) in the image
%domain from the C matrix in the k-space due to the iFFT
%   It receives as parameters:
%       - CK: the original C matrix in the k-space domain (stationary)
%       - Size: the size of the 2D-iFFT
%       - L: number of coils
%       - Lines: the lines acquired in the k-space
%
%   Output parameter:
%       - CX: C covariance matrix for every pixel in the image across
%       different channels.

%We compute the iFFT matrix
F=ifft(fftshift(eye(Size),1),[],1);

%We commpute the C matrix for transforming a line. Lines are fully sampled.
ownCorrRow=diag(F*F.');
nonZCol=find(abs(ownCorrRow)>max(abs(ownCorrRow))*10^-3);
valCol=ownCorrRow(nonZCol);
Z=zeros(Size,1);
Z(Lines)=1;
Rows=diag(Z);

%Then we compute the C matrix for the columns. This just applies to the
%columns that were non-zero in the previous step.
ownCorrCol=diag(F*Rows*F.');
nonZRow=find(abs(ownCorrCol)>max(abs(ownCorrCol))*10^-3);
valRow=ownCorrCol(nonZRow);

%We fill out the matrix containing the results
MX=zeros(Size,Size);
MX(nonZRow,nonZCol)=valRow*valCol.';
CX=repmat(permute(CK,[4 3 1 2]),[Size Size 1 1]).*repmat(MX,[1 1 L L]);
end

