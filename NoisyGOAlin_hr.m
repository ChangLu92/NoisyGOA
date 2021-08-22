function [] = NoisyGOAlin_hr(dataseth,datasetr,goObj,GOs)
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
fprintf('start %s  at %s\n,==Method:%s==',dataseth, datestr(now),'NoisyGOAlin_hr');
load(dataseth);
load(datasetr);

selGOs=GOs;
size_go=length(selGOs);
k=50;  % number of nearest neighbourhood proteins
% datapath=[pwd,filesep,'funsim20160409',filesep];%pwd is the current work directory
% addpath(datapath);
if size_go==3865
    gnd1=hGO.ccLabels;
    gnd2=rGO.ccLabels;
    rootGO=5575;   %ccroot
%     load('lin_sim_struct_cc.mat');
%     load('funsim_lin_struct_cc.mat');
end

if size_go==9982
    gnd1=hGO.mfLabels;
    gnd2=rGO.mfLabels;
    rootGO=3674;   %mfroot
%     load('funsim_lin_struct_mf.mat');
end

if size_go==28195
    gnd1=hGO.bpLabels;
    gnd2=rGO.bpLabels;
    rootGO=8150;  %bproot
%     load('funsim_lin_struct_bp.mat');
end
parGOs=getParentGOs(selGOs, selGOs,goObj);
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
lin_sim=zeros(Nfun,Nfun);
selGOProbs=sum(gnd_h,1)/Ndata;
% flag=sum(gnd3,2);
for ii=1:Nfun
%         fprintf('%d row compute linsim: %s\n', ii, datestr(now));
    go1=selGOs(ii);
    if(selGOProbs(ii)~=0)
    for jj=ii+1:Nfun
        go2=selGOs(jj);
        lin_sim(ii,jj)=compGOLinSim(go1,go2, parGOs{ii},parGOs{jj},selGOs, selGOProbs);
    end
    end
    lin_sim(ii,ii)=1;
    
end
lin_sim=max(lin_sim,lin_sim');
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

%%
r_score2=zeros(Ndata,Nfun);
for i=1:Ndata
    Idx=find(gnd_h(i,:)==0);
    knn_gnd=gnd_h(knn_ind(i,:),:);
    sum_knn_gnd=sum(knn_gnd);
    sum_knn_gnd(Idx)=0;
    sum_knn_gnd=sum_knn_gnd+gnd_h(i,:);
    b=sum_knn_gnd~=0;
    c=1./sum_knn_gnd(b);
    r_score2(i,b)=c;
end
[value2,ind2]=sort(r_score2,2,'descend');

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

% for ii=1:Ndata
%     if num_perprotein_noise(ii)>0
%     Idx=find(newgnd(ii,:)>0);  
% %     candidate=ind2(ii,1:num_perprotein_noise(ii)*1.5);
% %      candidate=ind2(ii,1:num_perprotein_noise(ii));
%      vv=value2(ii,num_perprotein_noise(ii));  %确定最后一个value的值
%      numnoise=length(find(value2(ii,:)>=vv));       %确定一共应该有多少candidate
%      candidate=ind2(ii,1:numnoise);
%     [~,ind3]=sort(r_score(ii,candidate),2,'descend');
% %     candidate=candidate(ind3);
%     noise=candidate(ind3(1:num_perprotein_noise(ii)));
% %     noise=ind(ii,1:num_perprotein_noise(ii));
%     childnoise=childGOs(noise);
%     newgnd(ii,noise)=0;
%     for jj=1:length(childnoise)
%         childidxp=getGOIdx(childnoise{jj},selGOs);
%         childIdx=intersect(childidxp,Idx);
%         newgnd(ii,childIdx)=0;
%     end
%     end
% end

%% compute precision, recall and f1-measure
Y=gnd_r;
Z =newgnd(1:Ndata,:);
[tp,per_pre,per_re,per_f1]=PRF(gnd_h,Y,Z);
data=length(find(num_perprotein_noise>0));
prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall，Maf1';
precision=cell(10,1);
precision{1}=tp;
precision{2}=per_pre;
precision{3}=per_re;
precision{4}=per_f1;
precision{5}=sum(per_pre)/data;
precision{6}=sum(per_re)/data;
precision{7}=sum(per_f1)/data;

if rootGO==5575
    evalstr=['save results',filesep,dataseth, '_NoisyGOAlin_hr_cc_k50.mat precision prec_seq'];
end
if rootGO==3674
    evalstr=['save results',filesep,dataseth, '_NoisyGOAlin_hr_mf_k50.mat precision prec_seq'];
end
if rootGO==8150
    evalstr=['save results',filesep,dataseth, '_NoisyGOAlin_hr_bp_k50.mat precision prec_seq'];
end
eval(evalstr);

fprintf('\n =====finish NoisyGOA_hr time=%s\n',datestr(now));




