%%%%%%%%%%%%%%%%%%%% MRMR����ѡ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
geneexpressionall = csvread('E:\LUAD_model\gene_id\gene_expression_all.csv'); %LUAD���еĻ�����ֵ���������ˣ�
d = geneexpressionall';
d=zscore(d);%����׼��
clusternortumor = csvread('E:\LUAD_model\gene_id\cluster.csv');%���ͱ�ǩ��0,1,2,3��
f = clusternortumor';
fea = mrmr_miq_d(d,f,500);%MRMR����ѡ��

%%%%%%%%%%%%%%%%%%%% IFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geneexpression = csvread('E:\LUAD_model\gene_id\gene_expression.csv');%LUAD��֢���ߵĻ�����ֵ�����������ˣ�
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