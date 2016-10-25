% Function: coilCombine
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,nz,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nz,1,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011

function [im2,myfilt] = coilCombine_SingleSlice(im1,myfilt)

% Let's make the coil dimension the fourth one and the TE the third
im1 = permute(im1,[1 2 4 3]);

% Get image dimensions and set filter size
[sx,sy,N,C] = size(im1);


% Initialize
im2 = zeros(sx,sy,N);
% Rs = zeros(sx,sy,C,C,N);

FlagFilt=isempty(myfilt);
if(FlagFilt)
    filtsize = 7;
    Rs = zeros(sx,sy,C,C);

    % Get correlation matrices
    for kc1=1:C
      for kc2=1:C
        for kn=1:N
    %       Rs(:,:,kc1,kc2,kn) = Rs(:,:,kc1,kc2,kn) + filter2(ones(filtsize),im1(:,:,kn,kc1).*conj(im1(:,:,kn,kc2)),'same'); %PARA IR ECO A ECO
          Rs(:,:,kc1,kc2) = Rs(:,:,kc1,kc2) + filter2(ones(filtsize),im1(:,:,kn,kc1).*conj(im1(:,:,kn,kc2)),'same'); %PARA HACER TODOS LOS ECOS JUNTOS
        end
      end
    end
    
    % Compute and apply filter at each voxel
    RsCell=permute(num2cell(permute(Rs,[3 4 1 2]),[1 2]),[3 4 1 2]);
    myfiltCell=cellfun(@(x) getFilterWalsh(x),RsCell,'UniformOutput',0);
    myfilt=cell2mat(myfiltCell);
else
    myfiltCell=num2cell(myfilt,[3]);
end

im1Cell=permute(num2cell(permute(im1,[3 4 1 2]),[1 2]),[3 4 1 2]);
im2Cell=cellfun(@(x,y) reshape((squeeze(x)')*(y.'),[1 1 N]),myfiltCell,im1Cell,'UniformOutput',0);
im2=cell2mat(im2Cell);

% In case the input data are single
if strcmp(class(im1),'single')
  im2 = single(im2);
end
end

function FiltWalsh=getFilterWalsh(Rs)
    [U,~] = svd(Rs);
    FiltWalsh = permute(-U(:,1),[3 2 1]);
end
