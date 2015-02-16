function [g,r] = rgp(x,y,k)
  f = aug(x,k);
  z = f(x);
  a = inv((eye(size(z,2))*1e-6+z'*z))*(z'*y);
  g = @(x0) f(x0)*a;
  r = g(x);
