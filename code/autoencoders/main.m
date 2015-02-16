% Load and split data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MNIST
d1 = 28;
d2 = 28;
d3 = 1;
load('data/mnist.mat')
y  = [double(train0);double(train1);double(train2);...
       double(train3);double(train4);double(train5);...
       double(train6);double(train7);double(train8);...
       double(train9)]./255;

% CIFAR
% IMPORTANT! You must generate cifar.mat as "cat cifar.mat.* > cifar.mat"
% d1 = 32;
% d2 = 32;
% d3 = 3;
% y=load('data/cifar.mat');
% y=double(y.d)./255;


%load('frey_rawface.mat')
%y = double(ff')./255;

p_tr  = 0.8;
p_va  = 0.1;
[n,D] = size(y);
n_tr  = floor(n*p_tr);
n_va  = floor(n*p_va);
n_te  = n-n_tr-n_va;
i     = randperm(n);
i_tr  = i(1:n_tr);
i_va  = i((n_tr+1):(n_tr+n_va+1));
i_te  = i((n_tr+n_va+2):end);
y_tr  = y(i_tr,:);
y_va  = y(i_va,:);
y_te  = y(i_te,:);
n     = n_tr;

% Construct model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d=[100]
  for k1=[2000]
    for k2=[2000]
      np = (k1*d+k2*size(y,2))/prod(size(y))*100;
      fprintf('TRAIN: %ix%i, %.1f%%. RP: %i, %i. MODEL SIZE: %.0f%%. LATENTS: %i... ',...
        n_tr,size(y,2),p_tr*100,k1,k2,np,d);
      tic;
      [g,x_tr] = rpca(y_tr,k1,d); 
      [f,r_tr] = rgp(x_tr,y_tr,k2);
      r_va     = f(g(y_va));
      err_va   = sum(sum((y_va-r_va).^2))/sum(sum((y_va).^2));
      fprintf('ERR: %f. TIME: %f.\n', err_va,toc);
    end
  end
end

r_te = f(g(y_te));

% Fit GMM and sample some digits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen;
