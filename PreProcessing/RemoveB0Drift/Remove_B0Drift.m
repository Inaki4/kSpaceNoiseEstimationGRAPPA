function []=Remove_B0Drift(filepath,filename)
%Remove_B0Drift Function to correct for the B0 drift in the data
%comparing the g-factors in the case when the ACS region is present for the
%real phantom
%   - filepath: path where the file containing the acquisition is storesd
%   - filename: file to preprocess and correct for B0drift
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

%We load the file
load([filepath filename]);

%We are going to use the point of k-space with maximum energy since it will
%have the highest SNR (assuming stationary noise) and it will be easier to
%estimate the phase.
maxValue=max(reshape(abs(K(:,:,1,1)),[],1));
[mxX,mxY]=find(abs(K(:,:,1,1)) == maxValue);


%We interpolate the phase of the data as a cubic function of time and
%remove it
y=squeeze(angle(K(mxX,mxY,1,:)));
x=(1:200)';
p = polyfit(x,y,3)
A=[x.^3 x.^2 x ones(size(x))];y2=A*p.';
ExpCorrection=repmat(exp(complex(0,reshape(-y2,[1 1 1 length(y2)]))),[128 128 8]);
K=K.*ExpCorrection;

%To ensure steady-state, we use just the last 100 realizations
K=K(:,:,:,end-99:end);

%We save the pre-processed data
save([filepath filename(1:end-4) '_B0DriftRemoved' filename(end-3:end)],'K');
end