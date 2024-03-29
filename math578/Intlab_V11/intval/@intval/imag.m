function c = imag(a)
%IMAG         Implements  imag(a)  for intervals (imaginary part)
%
%   c = imag(a)
%
%result zeros(size(a)) for real a,
%  otherwise real interval of imaginary part of a
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  faster check for rounding to nearest
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  c = a;

  if a.complex
    c.complex = 0;
    if isequal(a.rad,0)
      c.inf = imag(a.mid);
      c.sup = c.inf;
    else
      rndold = getround;
      setround(-1)
      c.inf = imag(a.mid) - a.rad;
      setround(1)
      c.sup = imag(a.mid) + a.rad;
      setround(rndold)
    end
    c.mid = [];
    c.rad = [];
  else
    if issparse(a.inf)
      [m n] = size(a.inf);
      c.inf = sparse([],[],[],m,n,0);
    else
      c.inf = zeros(size(a.inf));
    end
    c.sup = c.inf;
  end
  