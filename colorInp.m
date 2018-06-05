clear;
close all;
run vlfeat-0.9.20/toolbox/vl_setup;

%% read image
f  = imread('images/image_House256rgb','png');

% - resize original image
   s = 1000./min(size(f,1),size(f,2)); % resize scale
   if s<1
       f = imresize(f,s);
   end
fg = double(rgb2gray(f));
f  = double(f);
f = rgb2ycbcr(f);

[m,n,k] = size(f);
N       = m*n;

if k ~= 3
    disp('input is not color image');
end

%% sample data randomly
randIndex       = randperm(m*n);
sample_num      = floor(m*n*0.1); % 0.1 sample rate
id0             = randIndex(1:sample_num);
id1             = randIndex(sample_num+1:end);
id_matrix       = zeros(m,n);
id_matrix(id0)  = 1;

%% initialization
% -- generate the initial guess
fw = zeros(m*n,k);
for i = 1:k
    fw(id0+(i-1)*m*n) = f(id0+(i-1)*m*n);
    fw(id1+(i-1)*m*n) = mean(f(id1+(i-1)*m*n))+std(f(id0+(i-1)*m*n))*randn(size(id1));
end
fw = reshape(fw,m,n,k);
% -- parameters of patches
p1=5; % size of patch in x direction
p2=5; % size of patch in y direction

%% initailization
% -- initailize u as the intensity of I
id_p = image2patch(id_matrix,0,0,'center');

id   = find(id_p);
idi  = find(~id_p);

uf = zeros(m*n,k);
for i = 1:k
    uf(:,i) = image2patch(fw(:,:,i),0,0,'center');
end

u0 = image2patch(fg,0,0,'center');

x1=[1:1:m];
x2=[1:1:n];
[X,Y]=meshgrid(x1,x2);

X1=reshape(X,[],1);
X2=reshape(Y,[],1);

X1_c=X1*max(abs(u0(id)))/m;
X2_c=X2*max(abs(u0(id)))/n;

L1_c=3*max(fw(id))/m;
L2_c=3*max(fw(id))/n;

B = sparse(N,N);
D = sparse(N,N);

%% marching
% -- start looping
iterNum = 10;            % number of iterations
out     = cell(iterNum,1);   % output of every step result
lambda  = 10; 

disp('Iterating...');
tic

% v = VideoWriter('image_inpaint.avi','Uncompressed AVI');
% v.FrameRate=0.5;
% open(v);

for step = 1:iterNum 
    u0 = fw;
    % -- build patches from current image
    patches   = image2patch(u0, p1, p2, 'center');
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
    for i = 1:k
    
    u = uf(:,i);
    b = Ws*u(id) + WTs'*u(id) + sum(Q1.*SW1,2) - sum(Q2.*SW2,1)' - sum(Q3.*SW3,1)';
    L = ichol(coe_matrix);
    
    [u1]=pcg(coe_matrix,b,1e-6,100,L,L',u(idi));
    u(idi) = u1;
    u0 = reshape(u,n,m)';
    uu = u0';
    u = uu(:);
    uf(:,i) = u;
    fw(:,:,i) = u0;
    end
    %%
    % -- shrinkage
%     D_u  = sparse(id_row,id_col,swd.*(u(id_row)-u(id_col)),N,N);
%     DuB = D_u + B;
%     
%     DD = DuB(id_row+(id_col-1).*N);
%     D_next = zeros(1,N);
%     for j = 1:N
%         D_next(1,j) = norm(DD(:,j));
%     end
%     D_new = max(repmat(D_next,50,1)-lambda,0)./repmat(D_next,50,1);
% 
%     D = sparse(id_row,id_col,D_new,N,N).*DuB;
%     %%
%     % -- update the Lagrange multiplier
%     B = DuB - D; 
    
    fprintf('step=%d, PSNR=%f \n', step, psnr3d(ycbcr2rgb(fw),ycbcr2rgb(f)));
%     fprintf('step=%d, PSNR=%f \n', step, psnr3d(fw,f));
    out{step} = fw; % save the data
    imagesc(ycbcr2rgb(fw))
%     imagesc(uint8(fw))
%     un=(u0-min(min(u0)))/(max(max(u0))-min(min(u0)));
%     writeVideo(v,un);
    pause(0.2)
end
toc
disp(['Done in ' num2str(toc) ' secondes.']);
% close(v);
%% show results
figure; imshow(uint8(f)); axis image; 
figure; imshow(uint8(fg)); axis image; 
figure; imshow(uint8(fw)); axis image; 
figure; imshow((ycbcr2rgb(f/255))); axis image;
figure; imshow((ycbcr2rgb(fw/255))); axis image;
figure; imshow(id_matrix); axis image; 
figure;
for i = 1:10
    subplot(3,4,i); imshow(uint8(out{i}));
end