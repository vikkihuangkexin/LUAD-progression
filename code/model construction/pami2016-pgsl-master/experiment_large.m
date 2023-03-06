% filename = 'pendigits';
% 
% load(sprintf('%s.mat', filename));
X = csvread('E:\LUAD_model\gene_id\datacluster.csv',1,1);
X=X';
y=csvread('E:\LUAD_model\gene_id\cluster.csv');
ncls = length(unique(y));
ncenter = 100*ncls;

options.PCARatio = 0.95;
[P, ~] = PCA(X,options);
newX = X * P;

% normalization
max_dim = max(abs(newX),[],1);
newX = newX ./ repmat(max_dim, size(newX,1),1);

% k-means select centers
[N,D] = size(X);
rand('state',0);
[PI, centers]=litekmeans(newX,ncenter,'Replicates',10);

y_set = cell(ncenter,1);
for i=1:N
   y_set{PI(i)} = [y_set{PI(i)}, y(i)];
end

y_center = zeros(ncenter,1);
for i=1:length(y_center)
    tbl = tabulate(y_set{i});
    [max_val, max_idx] = max(tbl(:,2));
    y_center(i) = tbl(max_idx,1);
end

% run the regularized principal graph method
params = struct('maxiter', 20, ...
        'eps', 1e-5, ...
        'gstruct', 'span-tree',...
        'gamma', 0.1, ...
        'sigma', 0.1, ...
        'lambda', 0.1,...,
        'nn', 5,...
        'verbose',true);

C0=centers';
G =[];
if strcmp(params.gstruct,'l1-graph')
    C0 = centers';
    nC0 = size(C0, 2);
    if params.nn<nC0
        G = get_knn(C0, params.nn);
    else
        G = ones(nC0,nC0) - eye(nC0,nC0);
    end    
end

time = cputime;
[C, W, P,objs] = principal_graph(newX', C0, G, params);
fprintf('time cost=%f sec\n', cputime-time);

plot_3Dgraph(X', C, W, y_center); grid on;

% draw results
[val, idx] = sort(y_center,'ascend');
sort_y_center = y_center(idx);
W(W<1e-5) = 0;
sort_W = W(idx,idx);

figure1=figure; 
hold on;
spy(full(sort_W),20);
[~,idx,~] = unique(sort_y_center);
idx = [0; idx; length(sort_y_center)+1];
idx = idx(1:end)-1+0.5;

mz = meshgrid(idx,idx);
plot(idx,mz,'--k');
plot(mz,idx,'--k');

xlabel({''})