function [vpatches,patches] = image2patch(I, p1, p2, position)
% img example: I = imread('cameraman.jpg')
% p1, p2: patch size, =1*1, 3*3, 3*5
% patches: example: p1xp2xncolors matrix, ncolors: 3rd dimention of img
% vpatches: vectors of the patches. example: dxN matrix
% position: position of the patch 

%%
% - Default setting
if(~exist('position','var'))
    position = 'normal';
end

%% Input image
% - size of img
[nr, nc, ncolors] = size(I);
N = nr * nc;

%% Detect edge
% kappa = 30;
% dt = 0.01;
% im = I;
% c = anisodiffEdge(im, kappa, dt);
% u = c;
% u_d = (u-min(min(u)))/(max(max(u))-min(min(u)));

%% Build the patch data matrix
% disp('Building patch data matrix..');
% tic;

vpatches = zeros([N, p1*p2*ncolors]);
% for i=1:nr
%     for j=1:nc
%         vpatches(i+(j-1)*nr, p1*p2*ncolors+1) = 2.7*i;
%         vpatches(i+(j-1)*nr, p1*p2*ncolors+2) = 2.7*j;
%     end
% end

% imSz = size(I);
% xIdxs = 1:imSz(2);
% yIdxs = 1:imSz(1);
% patches = cell(length(yIdxs),length(xIdxs));
% midpatches = cell(length(yIdxs),length(xIdxs));
% trpatches = cell(length(yIdxs),length(xIdxs));

% - location of the patches
x1 = 1:1:nr;
x2 = 1:1:nc;
[xx1,xx2]=meshgrid(x1,x2);
X = reshape(xx1,[],1);
Y = reshape(xx2,[],1);



switch lower(position)
    case 'normal'
        % - right-down patches
        Iaugpost = zeros(nr+p1, nc+p2, ncolors);
        for i = 1:ncolors
            Iaugpost(:,:,i) = padarray(I(:,:,i), [p1,p2], 'symmetric', 'post');
        end
        for k=1:ncolors
            for i=1:p1
                for j=1:p2
                    vpatches(:,(k-1)*p1*p2+(i-1)*p2+j)=Iaugpost((k-1)*(nr+p1)*(nc+p2)+(Y+j-2)*(nr+p1)+(X+i-1));
                end
            end
        end
    case 'center'
        % - center patches
        Iaug = zeros(nr+2*p1, nc+2*p2, ncolors);
        for i = 1:ncolors
            Iaug(:,:,i) = padarray(I(:,:,i), [p1,p2], 'symmetric');
        end
        for k=1:ncolors
            for i=1:(2*p1+1)
                for j=1:(2*p2+1)
                    vpatches(:,(k-1)*(2*p1+1)*(2*p2+1)+(i-1)*(2*p2+1)+j)=Iaug((k-1)*(nr+2*p1)*(nc+2*p2)+(Y+j-2)*(nr+2*p1)+(X+i-1));
                end
            end
        end
    otherwise
        error('!invalid position')
end
% disp(['Done in ' num2str(toc) ' secondes.']);
% %%
% % - padding patched to full size
% for i = 1:length(yIdxs)
%     for j = 1:length(xIdxs)
%         for k = 1:ncolors
%             midpatches{i,j}(:,:,k) = padarray(patches{i,j}(:,:,k),[p1-size(patches{i,j},1), p2-size(patches{i,j},2)], 'post');
%         end
%         midpatches{i,j}(midpatches{i,j} == 0) = mean(mean(mean(patches{i,j})));
%     end
% end
% 
% %%
% % - reshape patches to vpatches
% for i = 1:length(yIdxs)
%     for j = 1:length(xIdxs)
%         for k = 1:ncolors
%             trpatches{i,j}(:,:,k) = midpatches{i,j}(:,:,k)';
%         end
%     end
% end
% for i = 1:length(yIdxs)
%     for j = 1:length(xIdxs)
%         for k = 1:ncolors
%             vpatches((i-1)*imSz(2)+j,:) = reshape(trpatches{i,j},[1, p1*p2*ncolors]);
%         end
%     end
% end
% weights = [1:floor(p1*p2*ncolors/2),ceil(p1*p2*ncolors/2):-1:1];
% weights = 100*p1*p2*ncolors:-100:1;
% weights = 1/sum(weights)*weights;
% weights = repmat(weights, N, 1);
% vpatches = vpatches.*weights;

% save vpatches.mat



% switch lower(position)
%     case 'normal'
%         % - right-down patches
%         for i = 1:length(yIdxs)
%             for j = 1:length(xIdxs)
%                 if u_d(i,j) < 0.1
%                     patches{i,j} = Iaug(yIdxs(i)+p1:yIdxs(i)+2*p1-1,xIdxs(j)+p2:xIdxs(j)+2*p2-1,:);
%                 elseif u_d(i,j) < 0.15
%                     patches{i,j} = Iaug(yIdxs(i)+p1:yIdxs(i)+2*p1-3,xIdxs(j)+p2:xIdxs(j)+2*p2-3,:);
%                 elseif u_d(i,j) < 0.2
%                     patches{i,j} = Iaug(yIdxs(i)+p1:yIdxs(i)+2*p1-5,xIdxs(j)+p2:xIdxs(j)+2*p2-5,:);
%                 elseif u_d(i,j) < 0.25
%                     patches{i,j} = Iaug(yIdxs(i)+p1:yIdxs(i)+2*p1-7,xIdxs(j)+p2:xIdxs(j)+2*p2-7,:);
%                 else
%                     patches{i,j} = Iaug(yIdxs(i)+p1:yIdxs(i)+2*p1-9,xIdxs(j)+p2:xIdxs(j)+2*p2-9,:);
%                 end
%             end
%         end
%     case 'center'
%         % - center patches
%         for i = 1:length(yIdxs)
%             Isub = Iaug(yIdxs(i)+p1-1-floor(p1/2):yIdxs(i)+p1-1+ceil(p1/2)-1,:,:);
%             for j = 1:length(xIdxs)
%                 patches{i,j} = Isub(:,xIdxs(j)+p2-1-floor(p2/2):xIdxs(j)+p2-1+ceil(p2/2)-1,:);
%             end
%         end
%     otherwise
%         error('!invalid position')
% end

