function [G, W] = get_knn(X, K)
% INPUT:
% X: DxN
% K: K neighbors of each data point excluding itself

N = size(X,2);
norm_sq = repmat(sum(X.^2),N,1);
dist_sq = norm_sq + norm_sq' - 2 .* X' * X;
[~, sort_idx] = sort(dist_sq, 2, 'ascend');
knn_idx = sort_idx(:,1:(K+1));

% an edge exists if two nodes are neighbors of each other
rows = reshape( repmat(knn_idx(:,1), 1, K), N*K, 1);
cols = reshape( knn_idx(:,2:(K+1)), N*K, 1);
Gtmp = sparse(rows, cols, ones(length(rows),1), N, N);
G = Gtmp + Gtmp';
G(G==2) = 1;
W = dist_sq .* G;
