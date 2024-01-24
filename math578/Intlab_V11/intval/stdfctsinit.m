function stdfctsinit
%STDFCTSINIT  Sets INTLAB constants for rigorous intval standard functions
%
%   stdfctsinit
%

% written  12/30/98     S.M. Rump
% modified 09/25/99     S.M. Rump  Various constants added, Matlab version
%                                  appended to stdfctsdata filename
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 08/21/04     S.M. Rump  file naming to cure Matlab/Unix bug
% modified 09/10/07     S.M. Rump  INTLAB path, file name
% modified 11/08/07     S.M. Rump  constant for log2 added
% modified 07/10/08     S.M. Rump  load path (thanks to Pete MacMillin, Colorado)
% modified 10/14/08     S.M. Rump  stdfctsfile path corrected
% modified 06/17/09     S.M. Rump  typo
% modified 08/26/12     S.M. Rump  INTLABPATH omitted
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
% modified 12/03/12     S.M. Rump  INTLABdata
% modified 05/31/13     S.M. Rump  data for error function
% modified 06/20/13     S.M. Rump  data for gamma function and modp
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 01/18/15     S.M. Rump  data for fft
% modified 10/11/15     S.M. Rump  data for psi function, global variables
% modified 01/16/16     S.M. Rump  improved error estimates for verifyfft
%

  global INTLAB_CONST

  setround(0)                           % set rounding to nearest

  load(intvalinit('INTLABdata'))
  INTLAB_CONST.INTVAL_POWER10 = INTLAB_INTVAL_POWER10;

  INTLAB_INTLABPATH = INTLAB_CONST.INTLABPATH;

  load([ INTLAB_INTLABPATH 'intval' filesep 'fft_data' ]);
  INTLAB_CONST.FFTDATA_D = d;
  INTLAB_CONST.FFTDATA_R = r;
  INTLAB_CONST.FFTDATA_NMAX = length(r);

  load([ INTLAB_INTLABPATH 'intval' filesep 'erfdata' ]);
  INTLAB_CONST.ERFCDATA = erfc_data;

  load([ INTLAB_INTLABPATH 'intval' filesep 'gammadata' ]);
  INTLAB_CONST.GAMMADATA = gamma_data;
  INTLAB_CONST.GAMMALNDATA = gammaln_data;
  INTLAB_CONST.GAMMAMIN = gamma_min;

  load([ INTLAB_INTLABPATH 'intval' filesep 'psidata' ]);
  INTLAB_CONST.PSIDATA = psi_data;

  % prime data for isregular with p = 3001171
  load([ INTLAB_INTLABPATH 'utility' filesep 'primedata' ]);
  INTLAB_CONST.PRIMEDATA = prime_data;

  %Assign precomputed data:

  %Fl-pt bounds for multiples of Pi and long representation for 1/Pi are computed:
  %
  % PI.PI       Bounds with
  %                 PI.PIINF <= pi <= PI.PISUP .
  % PI.TWOPI    Bounds with
  %                 TWOPIINF <= 2*pi <= TWOPISUP .
  % PI.PI2_4    Array p2_i of length 4 with correct digits for
  %                 pi/2 = sum( p2_i * beta^(-i) )        for i=1..4
  %             with  0 <= p2_i < beta = 2^25, beta < p2_1 < 2*beta
  % PI.PI2      Bounds with
  %                 PI.PI2INF <= pi/2 <= PI.PI2SUP .
  % PI.PI2      Bounds with
  %                 pi/2  in  PI.PI2MID +/- PI.PI2RAD
  % PI.PI4      Best lower bound with
  %                 PI.PI4INF <= pi/4 .
  % PI.PI2_k    Bounds for k=3 and k=4 with
  %                 sum(PI.PI2_i) + PI.PI2_(i+1)INF <= pi/2 <=
  %                 sum(PI.PI2_i) + PI.PI2_(i+1)SUP
  % PI.INV_PI2  Array phi_i of length digits+3 with correct digits for
  %                 2/pi = sum( phi_i * beta^(3-i) )        for i=1..digits+3
  %             with  0 <= phi_i < beta = 2^25
  % PI.REC2SQRT_PI bounds for 2/sqrt(pi)
  %

  INTLAB_STDFCTS_PI.PIINF = pow2(210828714,-26) + pow2(4467992,-51);
  INTLAB_STDFCTS_PI.PISUP = pow2(210828714,-26) + pow2(4467993,-51);
  INTLAB_STDFCTS_PI.TWOPIINF = pow2(421657428,-26) + pow2(4467992,-50);
  INTLAB_STDFCTS_PI.TWOPISUP = pow2(421657428,-26) + pow2(4467993,-50);
  INTLAB_STDFCTS_PI.PI2_4 = [52707178 17894214 2313292 13390192];
  INTLAB_STDFCTS_PI.PI2INF = pow2(105414357,-26) + pow2(4467992,-52);
  INTLAB_STDFCTS_PI.PI2SUP = pow2(105414357,-26) + pow2(4467993,-52);
  INTLAB_STDFCTS_PI.PI2MID = pow2(105414357,-26)+pow2(4467992,-52);
  INTLAB_STDFCTS_PI.PI2RAD = pow2(592202854,-83)+pow2(5337116,-108);
  INTLAB_STDFCTS_PI.PI4INF = pow2(52707178,-26) + pow2(71576856,-53);
  INTLAB_STDFCTS_PI.PI2_1 = 105414357*2^(-26);
  INTLAB_STDFCTS_PI.PI2_2 = 8935984*2^(-53);
  INTLAB_STDFCTS_PI.PI2_3INF = 74025356*2^(-80) + 51665926*2^(-106);
  INTLAB_STDFCTS_PI.PI2_3SUP = 74025356*2^(-80) + 51665927*2^(-106);
  INTLAB_STDFCTS_PI.PI2_3 = 74025356*2^(-80);
  INTLAB_STDFCTS_PI.PI2_4INF = 51665926*2^(-106) + 117912650*2^(-133);
  INTLAB_STDFCTS_PI.PI2_4SUP = 51665926*2^(-106) + 117912651*2^(-133);
  INTLAB_STDFCTS_PI.INV_PI2 = [ ...
           0           0           0    21361414    28915984    11096033 ...
     7699743    10918840     3594405    13409825    26231294    10667863 ...
    24833813    17676579    14918692    28990464    26364858    17098953 ...
     1900061    30816610    10919163    21238191    22211508     9477427 ...
    18614701     3129120    20550102     7538342    15192690     2292829 ...
    31429251    23330578    33552257    13371383    19861899    11801818 ...
     8238297    32930110    11574132    32032748    27760634    18266771 ...
    29018091    30139126    15998183    29119124    19834559    23852579 ...
    31676226     2818456     4652155 ];
  INTLAB_STDFCTS_PI.REC2SQRT_PIINF = hex2num('3ff20dd750429b6d');
  INTLAB_STDFCTS_PI.REC2SQRT_PISUP = hex2num('3ff20dd750429b6e');
  
  % inclusion of 2/sqrt(pi)
  INTLAB_STDFCTS_ERF.ONE_SQRTPIINF = hex2num('3fe20dd750429b6d');
  INTLAB_STDFCTS_ERF.ONE_SQRTPISUP = hex2num('3fe20dd750429b6e');
  INTLAB_STDFCTS_ERF.TWO_SQRTPIINF = hex2num('3ff20dd750429b6d');
  INTLAB_STDFCTS_ERF.TWO_SQRTPISUP = hex2num('3ff20dd750429b6e');


  % LOG2        The first 42 bits of log(2) such that
  %               LOG2.APP+LOG2.INF <= log(2) <= LOG2.APP+LOG2.SUP
  % LOG2_       Bounds for 1/log(2) with
  %                 LOG2_.INF <= 1/log(2) <= LOG2_.SUP
  % LOG10_      Bounds for 1/log(10) with
  %                 LOG10_.INF <= 1/log(10) <= LOG10_.SUP
  %

  INTLAB_STDFCTS_LOG2.APP = pow2(46516319,-26) + pow2(58530816,-52);
  INTLAB_STDFCTS_LOG2.INF = pow2(64908018,-70) + pow2(63399728,-97);
  INTLAB_STDFCTS_LOG2.SUP = pow2(64908018,-70) + pow2(63399729,-97);
  INTLAB_STDFCTS_LOG2_.INF = pow2(24204406,-24) + pow2(86737662,-52);
  INTLAB_STDFCTS_LOG2_.SUP = pow2(24204406,-24) + pow2(86737663,-52);
  INTLAB_STDFCTS_LOG10_.INF = pow2(29145009,-26) + pow2(86435086,-54);
  INTLAB_STDFCTS_LOG10_.SUP = pow2(29145009,-26) + pow2(86435087,-54);


  % E            It is  E.INF <= exp(1) <= E.SUP

  INTLAB_STDFCTS_E.INF = pow2(182420805,-26)+ pow2(18110313,-51);
  INTLAB_STDFCTS_E.SUP = pow2(182420805,-26)+ pow2(18110314,-51);
  
  % save data
  INTLAB_CONST.STDFCTS_EXP = INTLAB_STDFCTS_EXP;
  INTLAB_CONST.STDFCTS_E = INTLAB_STDFCTS_E;
  INTLAB_CONST.STDFCTS_LOG = INTLAB_STDFCTS_LOG;
  INTLAB_CONST.STDFCTS_LOG2 = INTLAB_STDFCTS_LOG2;
  INTLAB_CONST.STDFCTS_LOG2_ = INTLAB_STDFCTS_LOG2_;
  INTLAB_CONST.STDFCTS_LOG10_ = INTLAB_STDFCTS_LOG10_;
  INTLAB_CONST.STDFCTS_PI = INTLAB_STDFCTS_PI;
  INTLAB_CONST.STDFCTS_ERF = INTLAB_STDFCTS_ERF;
  INTLAB_CONST.STDFCTS_SIN = INTLAB_STDFCTS_SIN;
  INTLAB_CONST.STDFCTS_COS = INTLAB_STDFCTS_COS;
  INTLAB_CONST.STDFCTS_TAN = INTLAB_STDFCTS_TAN;
  INTLAB_CONST.STDFCTS_ATAN = INTLAB_STDFCTS_ATAN;
  INTLAB_CONST.STDFCTS_SINH = INTLAB_STDFCTS_SINH;
  INTLAB_CONST.STDFCTS_TANH = INTLAB_STDFCTS_TANH;
  INTLAB_CONST.STDFCTS_SUCCESS = INTLAB_STDFCTS_SUCCESS;
  