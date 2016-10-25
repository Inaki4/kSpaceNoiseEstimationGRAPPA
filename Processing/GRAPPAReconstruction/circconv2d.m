function [Y] = circconv2d(X,h)
%CIRCCONV2D It computes the circular convolution of X with filter h
%We get the size of the filter. Ideally both dimensions are odd so the
%filter can be centered in each point. If one dimension is even, we take
%the center in the first half, like if it is 4x4, we take the center at the
%point (2,2).
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

[hx,hy]=size(h);
cx=ceil(hx/2);cy=ceil(hy/2);
Y=conv2(X([(end-(cx-2)):end 1:end 1:(hx-cx)],[(end-(cy-2)):end 1:end 1:(hy-cy)]),h,'valid');
end

