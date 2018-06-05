clear;
close all;
run vlfeat-0.9.20/toolbox/vl_setup;

%% read img
I     = imread('barbara_256by256','png');
% I = rgb2gray(I);s
I = double(I);
[m,n] = size(I);
N = m*n; % data size
I_gt = I;

%% sample data randomly
randIndex       = randperm(m*n);
sample_num      = floor(m*n*0.01); % 0.1 sample rate
id0             = randIndex(1:sample_num);
id1             = randIndex(sample_num+1:end);
id_matrix       = zeros(m,n);
id_matrix(id0)  = 1;

%% initialization
% -- generate the initial guess
I_ini      = zeros(m*n,1);
I_ini(id0) = I(id0);
% fill in the missing pixels with random number
I_ini(id1) = mean(I_ini(id0))+std(I_ini(id0))*randn(size(id1)); 
I_ini      = reshape(I_ini,m,n);
% -- parameters of patches
p1=5; % size of patch in x direction
p2=5; % size of patch in y direction

%% initailization
% -- initailize u as the intensity of I
It = I_ini';
u = It(:);
id_matrixt = id_matrix';
id_p = id_matrixt(:);
id = find(id_p);
idi = find(~id_p);

x1=[1:1:m];
x2=[1:1:n];
[X,Y]=meshgrid(x1,x2);

X1=reshape(X,[],1);
X2=reshape(Y,[],1);

X1_c=X1*max(abs(u(id)))/m;
X2_c=X2*max(abs(u(id)))/n;

L1_c=3*max(I_gt(id))/m;
L2_c=3*max(I_gt(id))/n;

B = sparse(N,N);
D = sparse(N,N);

%% marching
% -- start looping
iterNum = 25;            % number of iterations
out     = cell(iterNum,1);   % output of every step result
lambda  = 1; 
disp('Iterating...');
tic
u0 = I_ini;
% v = VideoWriter('image_inpaint.avi','Uncompressed AVI');
% v.FrameRate=0.5;
% open(v);
% for step = 1:iterNum
err = 1;
error = zeros(1,25);
for step = 1:iterNum
    
    % -- build patches from current image
    uo = u0;
    patches   = image2patch(I_gt, p1, p2, 'center');
    data      = [patches,L1_c*X1_c,L2_c*X2_c];
    [id_row,id_col,w,~,sw,swd] = ann(data', 50, id); % W is weight matrix, SW is sqrt(W)
    %%
    % -- solve the linear system
    W  = sparse(id_row,id_col,w,N,N);
    SW = sparse(id_row,id_col,sw,N,N);
    r = m*n - sample_num;
    Wps = W(idi,idi);
    DW  = sparse(1:r,1:r,sum(Wps,2),r,r);
    DWT = sparse(1:r,1:r,sum(Wps,1)',r,r);
    Ws  = W(idi,id);
    WTs = 100*W(id,idi);
    DWs = sparse(1:r,1:r,sum(Ws,2),r,r);
    DWTs= sparse(1:r,1:r,sum(WTs,1)',r,r);
    coe_matrix = DW + DWT - Wps - Wps' + DWs + DWTs ;
    Q  = D - B;
    Q1 = Q(idi,:);
    SW1= SW(idi,:); 
    Q2 = Q(idi,idi);
    SW2= SW(idi,idi);
    Q3 = Q(id,idi);
    SW3= 10*SW(id,idi);
    b = Ws*u(id) + WTs'*u(id) + sum(Q1.*SW1,2) - sum(Q2.*SW2,1)' - sum(Q3.*SW3,1)';
    L = ichol(coe_matrix);
    [u1]=pcg(coe_matrix,b,1e-6,100,L,L',u(idi));
    u(idi) = u1;
    u0 = reshape(u,n,m)';
    uu = u0';
    u = uu(:);
    %%
    % -- shrinkage
    D_u  = sparse(id_row,id_col,swd.*(u(id_row)-u(id_col)),N,N);
    DuB = D_u + B;
    
    DD = DuB(id_row+(id_col-1).*N);
    D_next = zeros(1,N);
    for j = 1:N
        D_next(1,j) = norm(DD(:,j));
    end
    D_new = max(repmat(D_next,50,1)-lambda,0)./repmat(D_next,50,1);

    D = sparse(id_row,id_col,D_new,N,N).*DuB;
    %%
    % -- update the Lagrange multiplier
    B = DuB - D; 
    
    fprintf('step=%d, PSNR=%f \n', step, psnr(u0,I_gt));
    out{step} = u0; % save the data
    imagesc(u0)
    colormap('gray')
    err = norm(u0-uo)/sqrt(m*n)
    error(step) = err;
    un=(u0-min(min(u0)))/(max(max(u0))-min(min(u0)));
%     writeVideo(v,un);
    pause(0.2)
end
toc
disp(['Done in ' num2str(toc) ' secondes.']);
% close(v);
%% show results
figure; colormap gray;
subplot(221); imshow(I_gt,[]); axis image; title('Original image');
subplot(222); imshow(id_matrix,[]); axis image; title('Subseample image');
subplot(223); imshow(out{10},[]); axis image; 
title('Weighted nonlocal TV');
figure; colormap gray;
for i = 1:10
    subplot(3,4,i); imshow(out{i},[]);
end
figure; imshow(out{10},[]); axis image; 