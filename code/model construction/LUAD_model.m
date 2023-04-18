%%%%%%%%%%%%%%%%%%%%%%%% MRMR feature selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
geneexpressionall = csvread('E:\LUAD_model\gene_id\gene_expression_all.csv'); %LUAD gene expression matrix(with normal sample)
% geneexpression = csvread('E:\LUAD_model\gene_id\gene_expression.csv');%LUAD cancer gene expression(without normal)
d = geneexpressionall';
d=zscore(d);%normalization
clusternortumor = csvread('E:\LUAD_model\gene_id\cluster.csv');%subtype labels(0,1,2,3)
f = clusternortumor';
fea = mrmr_miq_d(d,f,314);%MRMR feature selection

% cluster10 = csvread('E:\LUAD_model\gene_id\kmeans_cluster.csv');%kmeans_cluster is result for kmeans
cluster101 = csvread('E:\LUAD_model\kmeans_cluster1.csv');
cluster10 = clusternortumor;
cluster4 = find(cluster101==5);
cluster10([cluster4],:)=[];%delete label4 sample
geneexpressionall(:,[cluster4])=[];

%%%%%%%%%%%%%%%%%%%%%%%  PCA visualization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = geneexpressionall(fea,:);
data = data';
data=zscore(data);
% writematrix(data,'E:\LUAD_model\gene_id\data3141.csv');
[COEFF,SCORE,latent,tsquared,explained,mu]=pca(data);% PCA
dataPCA=SCORE(:,1:3);% pick top 3 components
X= dataPCA;
y = cluster10';% 10 gene
[N,D] = size(X);
%%%%%%%%%%%%%%%%%%%% set parameters for principal curve %%%%%%%%%%%%%%%%%%%%%
params = struct('maxiter',100, ...
        'eps', 1e-5, ...
        'gstruct', 'span-tree',...
        'gamma', 0.004, ...
        'sigma', 15, ...
        'lambda', 1,...
        'nn',5,...
        'verbose',true);

C0=X;
G =[];
if strcmp(params.gstruct,'l1-graph')
    C0 = X;
    nC0 = size(C0, 2);
    if params.nn<nC0
        G = get_knn(C0, params.nn);
    else
        G = ones(nC0,nC0) - eye(nC0,nC0);
    end    
end

time = cputime;
%%%%%%%%%%%%%%%%%%% calculate principal curve %%%%%%%%%%%%%%%%%%%%%%%%
[C, W, P,objs] = principal_graph(X', C0', G, params);%C:centers for principal graph,W: principal graph matrix,C0: the initialization of centroids
                                                       %P: cluster assignment matrix
%%%%%%%%%%%%%%%%%%% draw principal curve %%%%%%%%%%%%%%%%%%%%%%%%%
W(W <1e-5) = 0;
[iidx, jidx, val] = find(sparse(W));
figure;
hold on;                                                       

for i=1:length(iidx)
    
    plot3( [C(1, iidx(i)), C(1, jidx(i))], [C(2, iidx(i)), C(2, jidx(i))],...
        [C(3, iidx(i)), C(3, jidx(i))], '-k','LineWidth',5); %5
    hold on;

end

%%%%%%%%%%%%%%%%%%%%%% draw data sample %%%%%%%%%%%%%%%%%%%%%%%
X=X';

cluster101([cluster4],:)=[];
cluster9 = find(cluster101==9);
cluster101([cluster9],:)=8;
y = cluster101;
cls = unique(y);
ncls = length(cls);%ÀàÊý

colors = distinguishable_colors(10);
%colors = ["#268AFF","#DC2B14","#5EB63A","#AF7575","#F5E866","#BCD693","#AFD7DB","#3D9CA8","#C34C17"];
names={};
colors(4,1) = 27/255;
colors(4,2) = 166/255;
colors(4,3) = 140/255;
mark = ["o","^","s","p"];
size = [9,6,10,11];
for c = 1:ncls
    idx = find(y==cls(c));
    names{c}=int2str(cls(c));
    for m = 1:length(idx)
        h(c) = plot3(X(1,idx(m)), X(2,idx(m)), X(3,idx(m)),mark(cluster10(idx(m))+1),'Color',...
            colors(c,:), 'MarkerSize', size(cluster10(idx(m))+1),'MarkerFaceColor',colors(c,:));% 5
        hold on;
    end
end
fprintf('time cost=%f sec\n', cputime-time);
%%%%%%%%%%%%%%%%%%%%%  project data sample onto curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist = zeros(1,N);
sortid = zeros(1,N);
for m=1:N
    for n=1:N
        dist(1,n)=norm(X(:,m)-C(:,n));
    end
    sortid(1,m) = find(dist == min(dist));
    
end
cluster101(find(cluster101==6),:)=5;
cluster101(find(cluster101==7),:)=6;
cluster101(find(cluster101==8),:)=7;
for c = 1:N
     plot3( [C(1, sortid(c)), X(1, c)], [C(2, sortid(c)), X(2, c)],...
        [C(3, sortid(c)), X(3, c)], 'Color',colors(cluster101(c)+1,:),'LineWidth',0.1);
     hold on;  
end
