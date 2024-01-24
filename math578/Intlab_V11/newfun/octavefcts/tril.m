function res = tril(A,k)
%TRIL         Functionality as tril with check for parameter outside matrix
%

% written  08/04/14  S.M. Rump
%

  global INTLAB_CONST

  if INTLAB_CONST.OCTAVE    
    if nargin==1
      k = 0;
    end    
    [m n] = size(A);
    k = min(max(k,-m),n);
  end
  
  res = builtin('tril',A,k);
