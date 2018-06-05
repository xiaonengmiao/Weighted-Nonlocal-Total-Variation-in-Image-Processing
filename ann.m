function [id_row,id_col,w,wd,sw,swd] = ann(data, number_s, id)
% data: dxN matrix, d: dimension of data, N: number of data
% number_s: number of neighbors 
% get a sparse graph Laplacian 
%% 
datasize = size(data);
N = datasize(2);

kdtree = vl_kdtreebuild(data);
[idx, dist] = vl_kdtreequery(kdtree, data, data, 'NumNeighbors', number_s, 'MaxComparisons', min(N,2^10));

sigma = sparse(1:N,1:N,1./max(dist(21,:),1e-2),N,N);

id_row = repmat(1:N,number_s,1);
id_col = double(idx);
w = exp(-(dist*sigma).^2);
sw = sqrt(exp(-(dist*sigma).^2));
swd = sw;
wd = sw;
rio = 10;
swd(:,id) = sw(:,id)*rio;
wd(:,id) = sw(:,id)*rio;

% %% W
% W=sparse(id_row,id_col,w,N,N);
% % W = max(W,W');
% % DW = sparse(1:N,1:N,sum(W,2),N,N);

