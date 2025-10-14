% Obsoleted problem formulation for Matlab ODE suite.

function [out1,out2,out3] = vdp1(t,y,flag)
if nargin < 3 | isempty(flag)
  out1 = [y(1).*(1-y(2).^2)-y(2); y(1)];
else
  switch(flag)
    case 'init'
      % Return tspan, y0, and options:
      out1 = [0 20];
      out2 = [2; 0];
      out3 = [];
    otherwise
      error(['Unknown request ''' flag '''.']);
  end
end
