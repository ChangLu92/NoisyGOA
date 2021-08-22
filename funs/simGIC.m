function [knsor,knind] = simGIC(Labels,selGOProbs,k)

[Ndata,Nfun] = size(Labels);
sim_gic=zeros(Ndata,Ndata);
sumUnion=zeros(Ndata,Ndata);
% selGOProbs=sum(Labels,1)/Ndata;
%sumMat =sum(Labels,2);
idx=find(selGOProbs~=0);
IC=-log2(selGOProbs(idx));
selGOProbs(idx)=IC;
MatIC=repmat(selGOProbs,Ndata,1);
MatA=Labels.*MatIC;
Intersection=MatA*Labels';
for ii=1:Ndata
%        fprintf('%d row compute simGIC: %s\n', ii, datestr(now));
       MatB=repmat(Labels(ii,:),Ndata,1);
%        Intersection=MatA.* Labels;
       Union=MatB+Labels;
%        sumInter=sum(Intersection.*MatIC,2);
%        sumInter=Intersection(ii);
       Union(find( Union>0 ))=1;
       sumUnion(:,ii)=sum(Union.*MatIC,2);
%        sim_gic(ii,:)=(sumInter./sumUnion)';
end
sim_gic=(Intersection./sumUnion)';
%sim_gic=sim_gic';%xt(i) to all the labeled samples is one row
[sorted,index]=sort(sim_gic,2,'descend');%对行向量排序
knsor=sorted(:,2:k+1);
knind=index(:,2:k+1);
fprintf('KNN is over !\n');