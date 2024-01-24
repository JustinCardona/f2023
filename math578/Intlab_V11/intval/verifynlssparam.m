function [List,ListS,ListData] = verifynlssparam(f,P,x0,opt,varargin)
%VERIFYNLSSPARAM  Nonlinear system parameter estimation
%
%For a function f:R^n->IR^m, find all x within the box x0 with f(x) in P. 
%Typically, this function is used for parameter estimation. Standard calls are
%
%   [ X , XS ] = verifynlssparam(f,P,x0)                or
%   [ X , XS , Data ] = verifynlssparam(f,P,x0)
%
%or with optional parameters (if opt is omitted or empty, default
%parameters are used)
%
%   [ X , XS ] = verifynlssparam(f,P,x0,opt)            or
%   [ X , XS , Data ] = verifynlssparam(f,P,x0,opt) 
%
%Then, on output, 
%   X   is an nxK interval matrix such that f(x) is included in P for all x in X(:,k)
%   XS  is an nxL interval matrix of candidates with that property. 
%For all x in x0 with f(x) in P, x is in one of the boxes X_i or XS_j.
%
%The output parameter Data stores data to continue the search, for example by
%
%   [ X , XS , Data] = verifynlssparam(f,P,x0)
%   for k=1:kmax
%     [ X , XS , Data ] = verifynlssparam(Data)
%   end
%   
%Similarly,
%
%   [ X , XS ] = verifynlssparam(f,P,x0,[],param1,param2,...)   or
%   [ X , XS , Data ] = verifynlssparam(Data,opt)
%
%evaluates the function f with extra parameters.
%
%Upon completion, X is an n x K array, so that f(x) is in P for all x in
%some column of X. Moreover, XS is an n x L array so that possibly f(x) is
%in P for x in some column of XS. For x neither in X nor in XS, surely f(x)
%is not in P.
%
%input    f         f:R^n->IR^m, function to find all x with f(x) in P
%         P         m-vector or interval vector
%         x0        box to be searched
%         opt.fields  optional, if not specified, default values are used
%             Display   0    no extra display (default)
%                       1    for 1<=n<=3, plots of the solution set
%                       inf  plot also discarded boxes
%             Boxes     64   each step bisection into how many boxes
%             Method    Evaluation method
%                       'intval'   standard interval evaluation
%                       'midrule'  standard interval evaluation together with
%                                     midpoint rule
%                       'slope'    slope evaluation
%                       'affari'   affine arithmetic evaluation
%             iFunMax   maximal number of function evaluations, default 1e5
%             TolXAbs   Termination absolute tolerance solution boxes, default 1e-5
%             TolXRel   Termination relative tolerance solution boxes, default 1e-3
%output   X         n x K array of K inclusion boxes of x with f(x) in P
%         XS        n x L array of L possible inclusion boxes of x with f(x) in P
%         Data      Data to continue search
%
%Parameters in opt may be set by "verifynlssparamset":
%   opt = verifynlssparamset('boxes',256);
%or directly
%   [X,XS] = verifynlssparam(@Griewank,infsup(0,1),infsup(-600,600)*ones(2,1), ...
%       verifynlssparamset('display','~'));
%   size(X), size(XS)
%An exhaustive search for zero sets is performed by
%   [X,XS] = verifynlssparam(@Griewank,0,infsup(-600,600)*ones(2,1));
%The List XS of possible inclusion boxes might be long. To collect boxes
%with a significant part in common, use 
%   XS = collectList(XS);
%For long lists that might take a considerable amount of computing time,
%depending on the structure of XS. Therefore, collection in not included
%here. For details and examples, see collectList.
%The function f needs to be in "vectorized form", i.e. [f(x) f(y)] and f([x y]) 
%must be equal. To achieve that, replace operations "o" by ".o" and replace 
%direct access to components x(i) by x(i,:). 
%If f is a function handle or a character string, this may be achieved 
%by funvec(f); if f is an inline function, it is converted automatically. 
%
%The algorithm is based on
%
% S.M. Rump. Mathematically Rigorous Global Optimization in Floating-Point 
%   Arithmetic. to appear in Optimization Methods & Software, 2018.
%

% written  02/27/17  S.M. Rump
% modified 07/06/17  S.M. Rump  Final check for vectorization: f = @(x) g(g(x))
%                                 may fail (thanks to Thomas Wanner)
% modified 07/30/17  S.M. Rump  comments and check for method
% modified 07/30/17  S.M. Rump  maximal number of function evaluations
% modified 12/12/17  S.M. Rump  check call with correct data
% modified 12/25/17  S.M. Rump  filename for multiple call
% modified 03/18/18  S.M. Rump  comment
% modified 04/17/18  S.M. Rump  spell check
%

  global INTLAB_NLSS
  global INTLAB_CONST
  
  rndold = getround;
  if rndold
    setround(0)
  end
  
  % ignore input out of range; NaN ~ empty
  INTLAB_CONST.RealStdFctsExcptnIgnore = 1;
  % make sure that is reset
  dummy = onCleanup(@()restore);
  
  refine = isstruct(f);
  
  if refine                             % refinement call
    INTLAB_NLSS = f;
    if ~isequal(INTLAB_NLSS.Filename,mfilename)    % check for data
      error(['Stored data was generated by ' INTLAB_NLSS.Filename])
    end
    x0 = INTLAB_NLSS.x0;
    Pinf = INTLAB_NLSS.Pinf;
    Psup = INTLAB_NLSS.Psup;
    List = INTLAB_NLSS.List;
    v = INTLAB_NLSS.ListS;
    ListS = [];
    if nargin>1                         % additional options
      opt = verifynlssparamset(P);
    end
    
  else                                  % first call
    
    if ( ~isintval(x0) ) && ( ~isaffari(x0) )
      error('Constraint box must be intval or affari')
    end
    if size(x0,2)~=1
      error('Input box must be column vector')
    end
    
    INTLAB_NLSS = [];
    INTLAB_NLSS.x0 = x0;
    initconstants;
    
    if nargin<5
      param = {};
    else
      param = varargin;
    end
    INTLAB_NLSS.param = param;
    
    INTLAB_NLSS.F = funvec(f,mid(x0),[],param);  % vectorize f
    INTLAB_NLSS.M = checkfun(INTLAB_NLSS.F,x0,param);       % check f is vectorized
    INTLAB_NLSS.N = size(x0,1);                             % number of unknowns
    P = intval(P);
    Pinf = P.inf;
    Psup = P.sup;
    INTLAB_NLSS.Pinf = Pinf;
    INTLAB_NLSS.Psup = Psup;    
    
    v = x0;
    List = [];
    ListS = [];
    
  end
  
  if ( nargin<4 ) || isempty(opt)
    if refine
      opt = INTLAB_NLSS.opt;
    else
      opt = verifynlssparamset;
    end
  else
    opt = verifynlssparamset(opt);
  end
  INTLAB_NLSS.opt = opt;
  
  INTLAB_NLSS.SEE = opt.Display;
  INTLAB_NLSS.METHOD = opt.Method;
  INTLAB_NLSS.BOXES = opt.Boxes;
  INTLAB_NLSS.IFUNMAX = opt.iFunMax;
  INTLAB_NLSS.TolXAbs = opt.TolXAbs;
  INTLAB_NLSS.TolXRel = opt.TolXRel;
  
  ifun = 0;                     % initialize counter
  param = INTLAB_NLSS.param;
 
  if INTLAB_NLSS.SEE
    if INTLAB_NLSS.N>3
      warning('Graphical output up to 3 unknowns')
      INTLAB_NLSS.SEE = 0;
    else
      figure
      X = x0 + midrad(0,0.1*rad(x0));
      if length(x0)==1
        axis([inf(X) sup(X) -.1 .1])
      elseif length(x0)==2
        axis([inf(X(1)) sup(X(1)) inf(X(2)) sup(X(2))])
      else
        axis([inf(X(1)) sup(X(1)) inf(X(2)) sup(X(2)) inf(X(3)) sup(X(3))])
      end
      hold on
    end
    if isstruct(f) && INTLAB_NLSS.SEE && ( ~isempty(List) )
      plotintval(List,'r');
    end
  end
  
  if ( nargin>=4 ) && isstruct(x0)
    error('parameters already specified.')
  end
  
  while 1
    % halve only variables for constraint optimization
    [vinf,vsup] = halve(v.inf,v.sup,false);
    K = size(v,2);        % dimension and number of boxes
    while K<INTLAB_NLSS.BOXES
      Kold = K;
      [vinf,vsup] = halve(vinf,vsup,false);
      K = size(vinf,2);
      if Kold==K          % there might be only point-interval boundary boxes
        break
      end
    end
    v = intval(vinf,vsup,'infsup');
    switch INTLAB_NLSS.METHOD
      case 'intval'       % standard interval evaluation
        if isempty(param)
          y = feval(INTLAB_NLSS.F,v);
        else
          y = feval(INTLAB_NLSS.F,v,param{:});
        end
        ifun = ifun + size(v,2);
      case 'midrule'      % standard interval evaluation together with midpoint rule
        mv = mid(v);
        if isempty(param)
          y = feval(INTLAB_NLSS.F,intval(mv));
          yg = feval(INTLAB_NLSS.F,gradient(v,'matrixofvectors'));
        else
          y = feval(INTLAB_NLSS.F,intval(mv),param{:});
          yg = feval(INTLAB_NLSS.F,gradient(v,'matrixofvectors'),param{:});
        end
        ifun = ifun + 2*size(v,2);
        if size(v,2)==1
          y = intersect( yg.x , y + yg.dx*(v-mv) );
        else
          ygdx = permute(yg.dx,[3 2 1]);
          yy = y;
          for i=1:size(y,1)
            yy(i,:) = yy(i,:) + sum(ygdx(:,:,i).*(v-mv),1);
          end
          index = isnan(yy(:));
          if any(index(:))
            yy(index) = intval(-inf,inf,'infsup');
          end
          y = intersect(yg.x,yy);
        end
      case 'slope'          % slope arithmetic evaluation
        if isempty(param)
          y = feval(INTLAB_NLSS.F,slopeinit(v.mid,v));
        else
          y = feval(INTLAB_NLSS.F,slopeinit(v.mid,v),param{:});
        end
        ifun = ifun + size(v,2);
        y = y.r;
      case 'affari'       % affine arithmetic evaluation
        if isempty(param)
          y = feval(INTLAB_NLSS.F,affari(v));
        else
          y = feval(INTLAB_NLSS.F,affari(v),param{:});
        end
        ifun = ifun + size(v,2);
      otherwise
        error('invalid option for Method')
    end
    % check empty intersection
    index = any( bsxfun(@gt,Pinf,y.sup) | bsxfun(@lt,Psup,y.inf) , 1 ) | any(isnan(y),1);
    if isinf(INTLAB_NLSS.SEE)
      if any(index(:))
        plotintval(v(:,index),'w');
      end
    end
    % check fully verified
    index1 = all( bsxfun(@le,Pinf,y.inf) & bsxfun(@ge,Psup,y.sup) , 1 );
    if any(index1(:))
      List = [ List v(:,index1) ];
      if INTLAB_NLSS.SEE
        plotintval(v(:,index1),'r');
      end
    end
    index = index | index1;
    if any(index(:))
      v(:,index) = [];
    end
    % check absolute tolerance
    index = all( diam(v)<=INTLAB_NLSS.TolXAbs , 1 );
    if any(index(:))
      ListS = [ ListS v(:,index) ];
      if INTLAB_NLSS.SEE
        plotintval(v(:,index),'y');
      end
      v(:,index) = [];
    end
    % check relative tolerance
    index = all( rad(v)<=INTLAB_NLSS.TolXRel*mid(v) , 1 );
    if any(index(:))
      ListS = [ ListS v(:,index) ];
      if INTLAB_NLSS.SEE
        plotintval(v(:,index),'y');
      end
      v(:,index) = [];
    end
    
    % final check
    if isempty(v) || ( ifun>INTLAB_NLSS.IFUNMAX )
      break
    end
    
  end  % main loop
  
  ListS = [ ListS v];
  if INTLAB_NLSS.SEE && ( ~isempty(v) )
    plotintval(v,'y');
  end
  
  INTLAB_NLSS.List = List;
  INTLAB_NLSS.ListS = ListS;
  ListData = INTLAB_NLSS;
  ListData.Filename = mfilename;
  
  setround(rndold)
  
end  % verifynlssparam


function m = checkfun(f,x,param)
%CHECKFUN     check parallelization of function
%output: dimension of image space
%

% written  01/04/16     S.M. Rump
%

  if nargin<3
    param = [];
  end
  
  % randomize x
  x = infsup(mid(x)-rand(size(x)).*rad(x),mid(x)+rand(size(x)).*rad(x));
  
  mx = mid(x);
  xinf = x.inf;
  xsup = x.sup;
  
  % take care of infinite components
  index = isinf(x.inf);
  if any(index)
    xinf(index) = -1;
  end
  index = isinf(x.sup);
  if any(index)
    xsup(index) = 1;
  end
  index = isinf(xinf) & isinf(xsup);
  if any(index(:))
    xinf(index) = -1;
    xsup(index) = 1;
    mx(index) = 0;
  end
  x = [ infsup(xinf,mx) infsup(mx,xsup)];
  
  if isempty(param)
    y = f(x);
  else
    y = f(x,param{:});
  end
  m = size(y,1);
  y1 = intval([]);
  if isempty(param)
    for i=1:size(x,2)
      y1 = [y1 f(x(:,i))];
    end
  else
    for i=1:size(x,2)
      y1 = [y1 f(x(:,i),param{:})];
    end
  end
  if any(relerr(y,y1)>2e-15)
    error('The input function seems to be not parallelized. Please use "funvec" or ensure that f([x x])=[f(x) f(x)].')
  end
  
end  % function checkfun

function [xinf,xsup] = halve(xinf,xsup,constraint)
% Wg. Octave identische Kopie in verifyglobal 
% (kein Aufruf von private Funktionen aus private Funktionen)
  % Simple choice to bisect component with largest diameter; other choices
  % like maximal f'*rad(x) are possible, however, with mixed results:
  % sometimes significantly better, sometimes significantly worse.
  % Number of boxes is doubled
  
  global INTLAB_NLSS

  [n,K] = size(xinf);
  
  if constraint         % For constraint optimization, only variables are halved
    M = INTLAB_NLSS.M;
    xxinf = xinf(1:n-M,:);
    xxsup = xsup(1:n-M,:);
    lastrowinf = xinf(end-M+1:end,:);
    lastrowsup = xsup(end-M+1:end,:);
    xinf = xxinf;       % original variables
    xsup = xxsup;
  else
    xxinf = xinf;
    xxsup = xsup;
  end
  xxinf(isinf(xxinf)) = -realmax;
  xxsup(isinf(xxsup)) = realmax;
  mx = 0.5*xxinf + 0.5*xxsup;     % is finite
  rx = mx - xxinf;
  
  [dummy,index] = max(rx,[],1);
  if constraint
    Index = index + (0:K-1)*(n-M);
  else
    Index = index + (0:K-1)*n;
  end
  
  xiinf = xxinf(Index);
  xisup = xxsup(Index);
  midx = mx(Index);
  radx = rx(Index);
  
%   if constraint
%     diamxi = xisup-xiinf;
%   end
  % mx = mid(x(i)); take care of infinite intervals
  % do not take exact midpoint; do not use sign(xi.mid), it may be zero 
  magic = 0.0303195520061965;
  midx = midx + magic*(2*(xiinf>0)-1).*radx;
  narrow = ( radx <= 1e-12*max(abs(xiinf),abs(xisup)) );
  
  % take care of wide intervals
  wide = ( radx > 100 );    % => 10<sqrt(radx)<radx
  if any(wide)
    xwinf = xiinf(wide);
    xwsup = xisup(wide);
    midwx = midx(wide);
    if xwinf>0              % positive wide interval
      midwx = xwinf + sqrt(radx(wide)) + magic;
    elseif xwsup<0          % negative wide interval
      midwx = xwsup - sqrt(radx(wide)) - magic;
    else
      index = ( xwinf<-1 ) & ( xwsup>1 );  % e.g. [-2,realmax]
      if any(index)         % wide zero interval
        midwx(index) = magic;
      end
    end
    midx(wide) = midwx;
  end

  if all(narrow)    % do not bisect small diameter or point boxes
    return
  elseif any(narrow)
    % append new halved boxes but not narrow
    x1inf = xinf;
    x1sup = xsup;
    x1inf(Index) = midx;
    notnarrow = ~narrow;
    x1sup(Index(notnarrow)) = midx(notnarrow);
    xinf = [ xinf x1inf(:,notnarrow) ];
    xsup = [ x1sup xsup(:,notnarrow) ];
    if constraint
      xinf(n-M+1:n,:) = [ lastrowinf lastrowinf(:,notnarrow) ];
      xsup(n-M+1:n,:) = [ lastrowsup lastrowsup(:,notnarrow) ];
    end
  else
    % append new halved boxes
    x1inf = xinf;
    x1sup = xsup;
    x1inf(Index) = midx;
    x1sup(Index) = midx;
    xinf = [ xinf x1inf ];
    xsup = [ x1sup xsup ];
    if constraint
      xinf(n-M+1:n,:) = [ lastrowinf lastrowinf ];
      xsup(n-M+1:n,:) = [ lastrowsup lastrowsup ];
    end
  end
end  % function halve

