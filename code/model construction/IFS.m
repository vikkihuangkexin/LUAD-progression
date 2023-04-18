%%%%%%%%%%%%%%%%%%%% MRMR特征选择 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
geneexpressionall = csvread('E:\LUAD_model\gene_id\gene_expression_all.csv'); %LUAD所有的基因表达值（带正常人）
d = geneexpressionall';
d=zscore(d);%做标准化
clusternortumor = csvread('E:\LUAD_model\gene_id\cluster.csv');%亚型标签（0,1,2,3）
f = clusternortumor';
fea = mrmr_miq_d(d,f,500);%MRMR特征选择

%%%%%%%%%%%%%%%%%%%% IFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geneexpression = csvread('E:\LUAD_model\gene_id\gene_expression.csv');%LUAD癌症患者的基因表达值（不带正常人）
data = geneexpression(fea,:);
data = data';
% data=zscore(data);

Q=zeros(1,500);
for i=1:500
    X=data(:,1:i);
    count=0;
    for j=1:533
        [~,idx] = pdist2(X,X(j,:),'cityblock','Smallest',2);
        if(f(idx(2))==f(j))
            count=count+1;
        end
    end
    Q(1,i)=count/length(f);    
end
[m,index]=max(Q);
plot(Q,'*');