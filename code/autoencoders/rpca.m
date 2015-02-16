function [g,y] = rpca(x,k,d)
  f       = aug(x,k);
  z       = f(x);
  m       = mean(z);
  s       = std(z);
  s(s==0) = 1;
  z       = bsxfun(@rdivide,bsxfun(@minus,z,m),s);
  opts.k  = d;
  [~,~,a] = irlba(z,opts);
  y       = z*a;
  S       = diag(1./sqrt(var(y))); % whiten
  y       = y*S;
  g       = @(x0) bsxfun(@rdivide,bsxfun(@minus,f(x0),m),s)*a*S;
