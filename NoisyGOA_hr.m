function [] = NoisyGOA_hr(dataseth,datasetr,goObj,GOs)
% Identifying Noisy Gene Ontology Annotations using Taxonomic and Semantic
% Similarity on Archived GOA
% Chang Lu, Jun Wang, Guoxian Yu, College of Computer and Information
% Science, Southwest University. Contact gxyu@swu.edu.cn, lucifer1992@email.swu.edu.cn
%%%%%%%%%%
%   goObj: gene ontology read from GO file (obo format)
%   GOs:  GO term id for each column of the protein-function association matrix
%   dataseth: Name of the dataset archived on 2015-05-26.
%   datasetr: Name of the dataset archived on 2015-12-07.
%%%%%%%%%%
fprintf('start %s  at %s\n,==Method:%s==',dataseth, datestr(now),'NoisyGOA_hr');
load(dataseth);
load(datasetr);

selGOs=GOs;
size_go=length(selGOs);
k=5;  % number of nearest neighbourhood proteins

    gnd1=hGO.ccLabels;
    gnd2=rGO.ccLabels;
    rootGO=5575;   %ccroot
    load('lin_sim_struct_cc.mat');
% childGOs=getChildGOs(selGOs, selGOs,goObj);
% parGOs=getParentGOs(selGOs, selGOs,goObj);
%    lin_sim= getLinStructSim(selGOs,childGOs,parGOs, size_go);

childGOs=getChildGOs(selGOs, selGOs,goObj);
gnd_h=full(gnd1);
gnd_r=full(gnd2);

%only test on annotated proteins
index=find(sum(gnd_h,2)==0);
gnd_h(index,:)=[];
gnd_r(index,:)=[];
[Ndata,Nfun]=size(gnd_h);

num_perprotein_noise=zeros(Ndata,1); %the number of noisy annotations of each protein
gnd=gnd_r-gnd_h;

%%  count the number of noisy annotations
for ii=1:Ndata
    idx=find(gnd(ii,:)==-1);
    num_perprotein_noise(ii)=length(idx);
end

%% Identifying Noisy Gene Ontology Annotations
r_score=zeros(Ndata,Nfun);
[~,knn_ind]=KNN(gnd_h,k);
for i=1:Ndata
    idx=find(gnd_h(i,:)>0);
    knn_gnd=gnd_h(knn_ind(i,:),:);
    for j=1:length(idx)
        idx_j=idx(j);
        r_score(i,idx_j)= WV(knn_gnd,idx_j,lin_sim);
    end
end
[~,ind]=sort(r_score,2,'descend');

 newgnd=gnd_h;
for ii=1:Ndata
    Idx=find(newgnd(ii,:)>0);   
    noise=ind(ii,1:num_perprotein_noise(ii));
    childnoise=childGOs(noise);
    newgnd(ii,noise)=0;
    for jj=1:length(childnoise)
        childidxp=getGOIdx(childnoise{jj},selGOs);
        childIdx=intersect(childidxp,Idx);
        newgnd(ii,childIdx)=0;
    end
end

%% compute precision, recall and f1-score
Y=gnd_r;
Z =newgnd(1:Ndata,:);
[tp,per_pre,per_re,per_f1]=PRF(gnd_h,Y,Z);
data=length(find(num_perprotein_noise>0));
prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall£¬Maf1';
precision=cell(10,1);
precision{1}=tp;
precision{2}=per_pre;
precision{3}=per_re;
precision{4}=per_f1;
precision{5}=sum(per_pre)/data;
precision{6}=sum(per_re)/data;
precision{7}=sum(per_f1)/data;

evalstr=['save results',filesep,dataseth, '_NoisyGOA_hr_cc_k5.mat precision prec_seq'];

eval(evalstr);

fprintf('\n =====finish NoisyGOA_hr time=%s\n',datestr(now));




