function [CrossCov] = myCrossCov(v)
%MYCROSSCOV computes the covariance matrix between real and imaginary
%components of given vector
%   Input parameter:
%   - v: complex vector.
%   
%   Output parameter: 
%   - CrossCov: cross covariance between real and imaginar components of v.
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

N=length(v);
x=real(v);
y=imag(v);
xc = bsxfun(@minus,x,sum(x,1)/N);  % Remove mean
yc = bsxfun(@minus,y,sum(y,1)/N);  % Remove mean
CrossCov = (xc.' * yc) / (N-1);

end

