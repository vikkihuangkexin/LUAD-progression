clc;
clear all;
close all;
% geneexpressionall = csvread('E:\LUAD_model\otherdata\GSE31210\data.csv'); %LUAD所有的基因表达值（带正常人）
% % geneexpression = csvread('E:\LUAD_model\otherdata\GSE26939_subtype.csv');%LUAD癌症患者的基因表达值（不带正常人）
% d = geneexpressionall';
% d=zscore(d);%做标准化
% clusternortumor = csvread('E:\LUAD_model\otherdata\GSE31210\subtype.csv');%亚型标签（0,1,2,3）,N*1
% f = clusternortumor';
% fea = mrmr_miq_d(d,f,110);%MRMR特征选择

%cluster101 = csvread('E:\LUAD_model\otherdata\GSE31210\kmeans_cluster.csv');%kmeans_cluster为kmeans的聚类结果
cluster10 = csvread('E:\LUAD_model\otherdata\GSE31210\cluster_iCluster.csv');%kmeans_cluster为kmeans的聚类结果
% cluster4 = find(cluster10==4);
% cluster10([cluster4],:)=[];%去除标签4
% geneexpressionall(:,[cluster4])=[];%去除基因表达矩阵中标签为4的样本
%%%%%%%%%%%%%%%%%%%%%%%  PCA可视化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = geneexpressionall(fea,:);
data= csvread('E:\LUAD_model\datasel.csv');
% cluster10 = csvread('E:\LUAD_model\otherdata\validation\subtype.csv');
data = data';
data=zscore(data);

% geneexpression = geneexpressionall(:,1:501);
% data = geneexpression(fea,:);
% data = data';
% data=zscore(data);
% writematrix(data,'E:\LUAD_model\otherdata\validation\datacluster.csv');

[COEFF,SCORE,latent,tsquared,explained,mu]=pca(data);%进行PCA降维
dataPCA=SCORE(:,1:3);%选取PCA的前三个主成分
% dataPCA=dataPCA';
% dataPCA = mapminmax(dataPCA);
% dataPCA=dataPCA';
X= dataPCA;
y = cluster10';%10个基因同源组
[N,D] = size(X);
%%%%%%%%%%%%%%%%%%%% 设置主曲线参数 %%%%%%%%%%%%%%%%%%%%%
params = struct('maxiter',100, ...
        'eps', 1e-5, ...
        'gstruct', 'span-tree',...
        'gamma', 0.015, ...
        'sigma', 20, ...
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
%%%%%%%%%%%%%%%%%%% 计算主曲线 %%%%%%%%%%%%%%%%%%%%%%%%
[C, W, P,objs] = principal_graph(X', C0', G, params);%C:centers for principal graph,W: principal graph matrix,C0: the initialization of centroids
                                                       %P: cluster assignment matrix
%%%%%%%%%%%%%%%%%%% 画主曲线 %%%%%%%%%%%%%%%%%%%%%%%%%
W(W <1e-5) = 0;
[iidx, jidx, val] = find(sparse(W));
figure;
hold on;                                                       

for i=1:length(iidx)
    
    plot3( [C(1, iidx(i)), C(1, jidx(i))], [C(2, iidx(i)), C(2, jidx(i))],...
        [C(3, iidx(i)), C(3, jidx(i))], '-k','LineWidth',5);
    hold on;

end

%%%%%%%%%%%%%%%%%%%%%% 画数据点 %%%%%%%%%%%%%%%%%%%%%%%
X=X';
y = cluster101;
cls = unique(y);
ncls = length(cls);%类数

colors = distinguishable_colors(ncls);
colors(4,1) = 27/255;
colors(4,2) = 166/255;
colors(4,3) = 140/255;
mark = ["o","^","s","p"];
size = [9,10,10,11];
names={};
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

%%%%%%%%%%%%%%%%%%%%%  将数据点投影到主曲线上 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist = zeros(1,N);
sortid = zeros(1,N);
for m=1:N
    for n=1:N
        dist(1,n)=norm(X(:,m)-C(:,n));
    end
    sortid(1,m) = find(dist == min(dist));
end

for c = 1:N
     plot3( [C(1, sortid(c)), X(1, c)], [C(2, sortid(c)), X(2, c)],...
        [C(3, sortid(c)), X(3, c)], 'Color',colors(cluster101(c)+1,:),'LineWidth',0.2);
     hold on;  
end


% for c = 1:ncls
%     idx1 = find(y==cls(c));
%     for m=1:length(idx1)
%         plot3( [C(1, idx1(m)), X(1, idx1(m))], [C(2, idx1(m)), X(2, idx1(m))],...
%         [C(3, idx1(m)), X(3, idx1(m))], '-','color',colors(c,:),'LineWidth',0.2);
%         hold on;   
%     end
% end

