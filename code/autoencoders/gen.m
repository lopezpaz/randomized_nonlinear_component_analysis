n_plot = 15;
rows   = 2;
j      = randperm(size(y_te,1),rows*n_plot);
yplots = [];
rplots = [];
plots  = [];

for this_row=1:rows
    yplots = [];
    rplots = [];
    for i=1:n_plot
        y_i = reshape(y_te(j(i+n_plot*(this_row-1)),:),d1,d2,d3);
        r_i = reshape(r_te(j(i+n_plot*(this_row-1)),:),d1,d2,d3);
        
        for k=1:d3
            y_i(:,:,k) = rot90(y_i(:,:,k)',4);
            r_i(:,:,k) = rot90(r_i(:,:,k)',4);
        end
        
        yplots = [yplots y_i];
        rplots = [rplots r_i];
    end
    
    plots = [plots;yplots;rplots];
end


% silvermann's rule-of-thumb
%l = .1;
%S = eye(size(x_tr,2))*l;
%
%;rplots.*0+255];
%
%for i=1:5
%  zplots = [];
%  z = f(mvnrnd(x_tr(randperm(size(x_tr,1),n_plot),:),S));
%  ii = randperm(size(x_tr,1),1);
%  z = f(mvnrnd(repmat(x_tr(ii,:),n_plot,1),S));
%  for j=1:n_plot
%    zplots = [zplots rot90(reshape(z(j,:),d1,d2,d3)',4)];
%  end
%  plots  = [plots;zplots];
%end

imshow(plots)

