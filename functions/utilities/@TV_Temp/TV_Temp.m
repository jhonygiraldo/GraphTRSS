function  res = TV_Temp()

%res = TV_Temp()
%
% Implements a difference operator along the time dimensin for dynamic MRI
%
% Ricardo Otazo 2008

res.adjoint = 0;
res = class(res,'TV_Temp');

