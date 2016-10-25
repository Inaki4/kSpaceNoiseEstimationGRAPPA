function [F,I,N,mssLines,knlLines,mssR] = recongrappa_multik( varargin );

% [F,I,N ] = recongrappa(Wk, Sk, Y, [opt args]);
% 
% opt args:
%  'v', {0,1}              - verbosity off/on
%  'kernel', string        - grappa kernel size, x by y. e.g. '2x5', '4x1'
%  'N', (xy)-by-num_patterns-L matrix   - grappa reconstruction parameters
%  'acs', array of indices - index locations to use for ACS lines
%  'radius', r             - reconstruction radius.  When GRAPPA is used to 
%                            improve data coverage for coil sensitivity
%                            estimates, only the center of k-space needs to
%                            be reconstructed. r \in [ 0 .. 1 ]
%
% 
% This is an implementation of     
% "Generalized autocalibrating partially parallel acquisitions"
%   by M. A. Griswold, P. M. Jakob, et. al..
%  Mag. Reson. Med.   47(6):1202-1210.  Jun 2002

%
% Copyright 2004-2007  Scott Hoge  (shoge at ieee dot org)
%
% All rights reserved. Licensed according to the GNU General Public Licence
% (see http://www.gnu.org)
%

% uncombined images are generated for each coil in the array, by
% filling in k-space lines.
%
%
% AUTO-SMASH: estimate k-space lines via linear fit from neighboring lines
%
%  coil 1     o * o * - * o *
%  coil 2     o * o * - * o *
%  coil 3     o * o * - * o *
%  coil 4     o * o * - * o *
%
%  comp data  o * o * x * o *
%
% GRAPPA: data from all coils used to patch data in one coil
%
%  coil 1     o * o * - * o *
%  coil 2     o * o * - * o *
%  coil 3     o * o * - * o *
%  coil 4     o * o * x * o *
%
%  *: acquired line, o: unacquired line, -: auto-callibration line
%  x: reconstructed line
%
%  S_j (k_y - m\Delta k_y ) = \sum_{l=1}^L \sum_{b=0}^{N_b-1}
%                                  n(j,b,l,m) S_l(k_y - b A \Delta k_y)
%
% in a variable density sampling scheme, the center of k-space is fully
% sampled and can be used to determine the linear fit parameters.
%
% A   - acceleration factor
% N_b - number of blocks used in the reconstruction
% 
% N_b = 1 ----> GRAPPA = VD-AUTO-SMASH

Wk = varargin{1};
s  = varargin{2};
Y  = varargin{3};

acs = []; N = []; verbose=0; kernel = '2x3'; dks = []; prnt = 1;

vararg = varargin(4:end);
for cnt=1:2:length(vararg),
  if strcmp( vararg{cnt}, 'v' ),
    verbose = cell2mat(vararg(cnt+1));
  elseif strcmp( vararg{cnt}, 'N' ),
    N = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'dks' ),
    dks = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'kernel' ),
    kernel = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'acs' ),
    acs = vararg{cnt+1};
  elseif (strcmp( vararg{cnt}, 'report' ) | strcmp( vararg{cnt}, 'r' )),
    prnt = vararg{cnt+1};
  end;
end;

if length(size(Wk)) == 3, 
  [Wm,Wn,Wc] = size(Wk);
else,
  Wm = Wk(1); Wn = Wk(2); Wc = Wk(3);
end;

%%% correct for the toggling in acquisition...
% s = tmult( s, diag( (-1).^Y ), 1 );

%% zero-pad the k-space data
Sk = zeros( Wm, Wn, Wc );
Sk(Y,:,:) = s;

  if verbose, fprintf(['* ' mfilename ': kernel size = ' kernel '\n']); end;

d = diff(Y(:));
if ( isempty(acs) ),
  prfl = max( sum(abs( Sk(:,:,:)),3).' );

  if length( d == 1) == 1,
    acs = find( prfl == max(prfl(Y(find(d==1)))) );
  else
    acs = Y([ min(find(d==1)); find( d == 1)+1 ]);
  end;
  if verbose
      fprintf(['* ' mfilename ': acs =']);
      fprintf(' %d', acs );
      fprintf('\n');
  end
end;

if ( isempty(dks) ),
  %% determine which delta-k's to repair.
  dks = [];
  for cnt=max(d):-1:2,
    if (sum( cnt == d ) > 0),
      dks = [ cnt; dks ];
    end;
  end;
end;

% set the kernel spacings based on 'kernel' string
nky = str2double(kernel(1));               % should be even
nkx = str2double(kernel(3));               % should be odd

indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];

pattern = {};

%% foreach dk spacing %Each line can be reconstructed in several ways using a slide window approach (In ky-1 ways)
for cnt=1:length(dks),
  tmp = [ repmat([ '*' repmat('o',1,dks(cnt)-1) ],1,nky-1) '*' ];
  z2 = strfind(tmp,'o');
  [tmp2,indx] = sort( abs(z2 - z2( floor(median(1:length(z2))) )) );
  for cnt2 = 1:length(z2), % (dks(cnt)-1), % 
    cnt3 = indx(cnt2);
    pattern{length(pattern)+1} = tmp;
    pattern{length(pattern)}(z2(cnt3)) = 'x';
  end;
end;

% pattern{1} = '*xoo*o*';
% pattern{2} = '*oxo*o*';
% pattern{3} = '*oox*o*';

% find the recon patterns we can use 
plist = cell(length(pattern),1);

for cnt = 1:length(pattern),
  v = extrapolate_pattern(pattern{cnt});     

  % if cnt == 14, keyboard; end;

  for cnt2 = 1:length(acs);
    %% check if the GRAPPA pattern...
    indx = acs(cnt2) + v;

    %% can be accomodated by the sampling pattern 
    tmp = find( Y(:)*ones(1,length(indx)) == ones(length(Y),1)*indx );
    if ( length(tmp) == length(indx) ),
      plist{cnt} =  [ plist{cnt} acs(cnt2) ];
    end;
  end;
end;

if verbose, disp(['* ' mfilename ': sampling patterns in use = ' ]); end;
if verbose
    for cnt1=1:length(pattern);
        fprintf('  %s: %d acs lines =', pattern{cnt1}, length(plist{cnt1}) );
        fprintf(' %3d', plist{cnt1} );
        fprintf('\n');
    end
  % for cnt2 = 1:length( plist{cnt1} ),
  %   fprintf('  %s', pattern{plist{cnt1}} );
  %   % v = extrapolate_pattern(pattern{plist{cnt1}(cnt2)});
  %   % norm(vec(Sk( acs(cnt1)+v, :, : )));
  % end;
end;

if (length(plist) == 0), 
  error(['sorry, acquisition pattern is incompatible with GRAPPA recon' ...
         ' patterns']);
end;

st = clock; % tic;

if verbose, fprintf('\n'); end;

mss = zeros( size(Sk,1), 1 );         % set the missing lines index vector 
mss(Y) = 1;

if isempty(N),
  % N = zeros( length(indx3)*nky*size(s,3),length(pattern),size(s,3));

  if verbose, fprintf('calc recon param\n'); end;
  for cnt = 1:length(pattern),
    %% and for each pattern ... 
    %% A2 = [];     b2 = [];
    A = zeros( length(plist{cnt})*size(s,2), length(indx3)*nky*size(s,3) ); %This matrix has storage space for [(Nacslines-2)·length(line),NelementsKernel·Ncoils], so probably it's gonna keep the estimated kernels for each coil for each pixel in the ACS lines. We don't consider first and last ACS lines because they don't have the inmediate above or under line to be estimated from.

    if verbose, fprintf('%3d  %s: ', cnt, pattern{cnt}); end;

    v = extrapolate_pattern(pattern{cnt});
    N(length(N)+1).pattern = pattern{cnt};

    n = zeros(size(A,2),size(s,3));

    for cnt2 = 1:length(plist{cnt}),
      %% and ACS line ...

      %% if the sampling pattern accomodates the GRAPPA pattern...
      indx = plist{cnt}(cnt2) + v;

      %% determine the k-space filling coefficients 
      for cnt3=1:length(indx3),
        y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
        y0( y0 == 0 ) = size(Sk,2);
        
        A( size(s,2)*(cnt2-1) + (1:size(s,2)), ...
           cnt3:length(indx3):size(A,2) ) = ...
            reshape( permute(squeeze( Sk(indx,y0,:) ),[ 2 3 1 ]), ...
                     size(Sk,2), size(Sk,3)*length(indx) ); %We create a matrix where with kenel(x) groups of L columns, in which each column in a group is one of the lines for each coil. Moving to another group means moving to another line involved in the kernel. This lines are displaced in the phase encoding direction depending on k(y). So in the final matrix A, the groups are oo kernel(y) columns
      end;
      %% for an N-by-1 kernel, these next lines are enough...
      % A2 = [ A2; ...
      %         reshape( permute(squeeze( Sk(indx,:,:) ),[ 2 3 1 ]), ...
      %                  size(s,2), size(s,3)*length(indx) )
      %       ];
      % 
      % b2 = [b2; squeeze( Sk( plist{cnt}(cnt2), :, l ) ).' ];
    end;
    
    % Apinv = pinv(A);
    AA = A'*A;
    for l=1:size(Sk,3),
      if verbose, fprintf('.'); end;
      % A is the same for each coil ...
      b = vec( squeeze(Sk( plist{cnt}, :, l )).' );
      
      % n(:,l) = A \ b;
      % n(:,l) = Apinv * b;
      n(:,l) = cgsolv( AA, A'*b, zeros(size(A,2),1), size(A,2) );

      if (verbose>2), pltcmplx( A*n(:,l), b ); keyboard; end;
    end;
    N(length(N)).n = n;
    
    if verbose, fprintf('\n'); end;
  end;  
end;

if verbose, fprintf('recon missing lines\n'); end;
%% now go recon the missing lines of k-space:
%This loop can be improved by going only through the lines that are missing
%instead of checking all of them and seeing if they are acquired or not.
mssi=0; %============================================================ADDED===============================================
mssLines=cell(sum(1-mss),1); %============================================================ADDED=========================
knlLines=zeros(sum(1-mss),1); %============================================================ADDED=========================
mssR=zeros(sum(1-mss),1); %============================================================ADDED=========================
for cnt2=find(mss==0)',
    mssi=mssi+1; %============================================================ADDED======================================
  if verbose, fprintf(' line %3d, from lines ',cnt2); end;

  tmp = [];
  %% determine which recon pattern we can use:
  for cnt4 = 1:length(pattern),
    v = extrapolate_pattern(pattern{cnt4});
    indx = cnt2 + v;
    %     tmp = find( Y(:)*ones(1,length(indx)) == ones(length(Y),1)*indx );
    %     if  ( ... % ( (min(indx) > 1) & (max(indx) < size(Sk,1)) ) & ...
    %         ( length(tmp) == length(indx) ) ),
    %       break;                          % jump out of the loop
    %     end;
    
    %If we are at the first or last line, we can't estimate them because
    %there is no line beneath or above them.
    if ( indx(1) < 1) | (indx(end) > size(Sk,1)); continue; end;

    tmp = (mss( indx ) == ones(length(indx),1));
    if ( sum(tmp) == length(tmp)  ),
      break; 
    end;

  end;  
  
  if ( ( (min(indx) < 1) | (max(indx) > size(Sk,1)) ) | ...
       ( sum(tmp) ~= length(indx) ) ),
    
    if verbose
        fprintf('o');
        fprintf('\n'); 
    end
    continue; 
  end;
 
  
  %% recon k-space line based on that pattern:  
  A = zeros( size(s,2), length(indx3)*length(v)*size(s,3) );
  for cnt3=1:length(indx3),
    y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
    y0( y0 == 0 ) = size(Sk,2);
    
    A( (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
        reshape( permute(squeeze( Sk(indx,y0,:) ),[ 2 3 1 ]), ...
                 size(Sk,2), size(Sk,3)*length(indx) );
  end;
    
  % A2 = reshape( permute(squeeze( Sk(indx,:,:) ),[ 2 3 1 ]), ...
  %              size(Sk,2), size(Sk,3)*length(indx) );
  
  if (verbose), % ( sum(Sk(Y(mss(cnt2))+1,:,l)) ~= 0 ); 
    if verbose
        fprintf(' %3d',[ indx(:) ]);
        fprintf(' (pattern: %s)', pattern{cnt4} );
        fprintf('\n');
    end
    if (verbose>1),
      figure(1); imagesc( log10(abs(Sk(:,:,l)) + 1e-1) ); pause(0.5);
      keyboard; 
    end;
  else,
    if verbose, 
        fprintf('x');
    end
  end;

%   for cnt3=1:length(N),
%     if strcmp( N(cnt3).pattern, pattern{cnt4} ),
%       break;                            % jump out of the loop
%     end;
%   end;
  n = N(cnt4).n;
  %========================================================================
  %========================================================================
  %=================================ADDED==================================
  %========================================================================
  %========================================================================
  %We save the lines in the kernel for the current missing line and the
  %kernel used. We also save the acceleration for the current line
  IndxPattern=strfind(pattern{cnt4},'*');
  mssR(mssi,1)=IndxPattern(2)-IndxPattern(1);
  mssLines{mssi}=indx;
  knlLines(mssi)=cnt4;
  %========================================================================
  %========================================================================
  %========================================================================
  %========================================================================
  
  for l=1:size(Sk,3),
    % for each coil ...
    
    % fprintf(' line %3d, pattern %3d\n', cnt2, plist(cnt3) );
    Sk( cnt2, :, l ) = A * n(:,l);
  end;

  if (verbose>1), keyboard; end;
end;
if verbose, fprintf('  \n'); end;


% Sk(1:2:size(Sk,1),:,:) = -Sk(1:2:size(Sk,1),:,:);
% Sk(:,1:2:size(Sk,2),:) = -Sk(:,1:2:size(Sk,2),:);

F = Sk;
I = sqrt( sum(abs(ifft2(Sk)).^2,3) ); 
I = I./(prod(size(I)));
mask = 0;

rpttime = etime(clock,st);
if verbose
    fprintf(' GRAPPA recon time: %g (avg per kernel: %g)\n',  ... 
        rpttime, rpttime/length(pattern) ); 
end

function v = extrapolate_pattern(p)
 
r = find( p == 'x' );
v = find( p == '*' ) - r;

