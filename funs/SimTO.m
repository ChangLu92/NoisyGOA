function [knsor,knind]=SimTO(labels,k)

epsilon=eps;
[Np,~]=size(labels);

sumMat=sum(labels,2);
repM=repmat(sumMat,1,Np);

tempMat=labels*labels';

minMat=min(repM,repM');
idx=minMat>epsilon;
minMat(idx)=1./minMat(idx);

TOsim=tempMat.*minMat;

[sorted,index]=sort(TOsim,2,'descend');%对行向量排序
knsor=sorted(:,2:k+1);
knind=index(:,2:k+1);
fprintf('KNN is over !\n');