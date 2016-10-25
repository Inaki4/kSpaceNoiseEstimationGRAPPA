function wk=coef_conv(coef,kernel,tasa)

% Convolution matrix for GRAPPA in k-domain
% 
%  INPUTS:
% 		coef: grappa_coefs
% 		tasa: acceleration rate
%%
%
% PARALLEL MRI TOOLBOX
%
% Santiago Aja-Fernandez, LPI
% www.lpi.tel.uva.es/~santi
% Valladolid, 28/05/2012
%------------------------------------------
if length(size(coef))==2
	coef=ajusta_coef(coef);
end


[X,Y,Z,W]=size(coef);

if Z~=W
	error('Wrong nunber of coeffs');
end


wk=zeros([kernel(1)*tasa-1,Y,Z,Z]);
Ind=reshape(repmat((1:tasa:(kernel(1)*tasa-1)),tasa-1,1),[],1)+repmat((0:(tasa-2))',kernel(1),1);
wk(Ind,:,:,:)=coef;
wk(kernel(1)*tasa/2,(Y+1)/2,:,:)=reshape(diag(ones(1,Z)),[1 1 Z Z]);

%We need to flip the kernel in both dimensions
wkCell=num2cell(wk,[1 2]);
wkCell=cellfun(@(x) fliplr(flipud(x)),wkCell,'UniformOutput',0);
wk=cell2mat(wkCell);


% for kk=1:Z
% 
% %Matriz de convolución
% M_conv=zeros([kernel(1)*tasa-1,Y,Z]);
% 
%     for jj=1:Z
%         
%         Ind=reshape(repmat((1:tasa:(kernel(1)*tasa-1)),tasa-1,1),[],1)+repmat((0:(tasa-2))',kernel(1),1);
%         M_conv(Ind,:,jj)=coef(:,:,jj,kk);
%         
% %         for gg=1:X/2
% %             for hh=1:tasa-1
% %                 M_conv(2*gg-1,:,jj)=coef(gg,:,jj,kk);
% %                 M_conv((end-(2*gg-2)),:,jj)=coef(end-(gg-1),:,jj,kk);
% %             end
% %         end
% 
%         M_conv(kernel(1)*tasa/2,(Y+1)/2,kk)=1;
% 
%         wk(:,:,jj,kk)=M_conv(:,:,jj);
%     end
% end
end



