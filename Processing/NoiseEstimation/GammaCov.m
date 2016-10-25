function Gamma = GammaCov(x)
%GAMMACOV computes the Γ covariance matrix of a complex normal
%distribution, given by Γ=E[(Z-μ)·(Z-μ)'] from the given vector
%   Input parameter:
%   - x: complex vector.
%   
%   Output parameter: 
%   - Gamma: Γ covariance matrix of v.
%
%   Copyright (C) 2016 Iñaki Rabanillo <irabvil@lpi.tel.uva.es>
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

[m,n] = size(x);
xc = bsxfun(@minus,x,sum(x,1)/m);  % Remove mean
Gamma = (xc.' * conj(xc)) / (m-1);

end