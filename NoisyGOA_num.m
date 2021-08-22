function [] = NoisyGOA1(dataset,goObj,GOs)
% Identifying noisy gene ontology annotations using taxonomic and semantic
% similarity by adding fixed number of annotations
% Chang Lu, Jun Wang, Guoxian Yu, College of Computer and Information
% Science, Southwest University. Contact gxyu@swu.edu.cn, lucifer1992@email.swu.edu.cn
%%%%%%%%%%
%   goObj: gene ontology read from GO file (obo format)
%   GOs:  GO term id for each column of the protein-function association matrix
%   dataset: Name of the dataset, 'YeastGOA' means Yeast proteins annotated with terms.
%%%%%%%%%%
fprintf('start %s  at %s\n,==Method:%s==',dataset, datestr(now),'NoisyGOA1');
load(dataset);
selGOs=GOs;
size_go=length(selGOs);

if size_go==3865
    gnd=hGO.ccLabels;
    rootGO=5575;   %ccroot
    load('lin_sim_struct_cc.mat', 'lin_sim');
end

if size_go==9982
    gnd=hGO.mfLabels;
    rootGO=3674;   %mfroot
    load('funsim_lin_struct_mf.mat', 'lin_sim');
end

if size_go==28195
    gnd=hGO.bpLabels;
    rootGO=8150;  %bproot
    load('funsim_lin_struct_bp.mat', 'lin_sim');
end

%% parameter setting
num_noise=1:5;    % number of noisy annotations
round=5;         % time of the independent repeated experiment
k=5;              % number of nearest neighbourhood proteins

%% only test on annotated proteins
gnd2=full(gnd);
index=find(sum(gnd2,2)==0);
gnd2(index,:)=[];
[Ndata,Nfun]=size(gnd2);

tp=zeros(length(num_noise),round,Ndata); %the number of correctly predicted noisy annotations of every protein
per_pre=zeros(length(num_noise),round,Ndata); % precision of identifying noisy annotations of every protein
per_re=zeros(length(num_noise),round,Ndata);  % recall of identifying noisy annotations of every protein
per_f1=zeros(length(num_noise),round,Ndata);  % f1-measure of identifying noisy annotations of every protein
Macropre=zeros(length(num_noise),round);    % average precision
Macrore=zeros(length(num_noise),round);     % average recall
Macrof1=zeros(length(num_noise),round);     % average f1-measure

childGOs=getChildGOs(selGOs, selGOs,goObj);
DirectchildGOs=getDirectChildGOs(selGOs, selGOs,goObj);

for num_noise=1:5
    for run=1:round
        %% add noisy annotations
        num_noise_Mat=repmat(num_noise,Ndata,1);
        [gnd3, num_perprotein_noise]=addNoise1(gnd2,selGOs,num_noise_Mat,DirectchildGOs);
        fprintf('==finish addNoise1!\n');
        
        %% Identifying Noisy Gene Ontology Annotations
        Y=gnd2;
        t0=clock;
        r_score=zeros(Ndata,Nfun);
        [~,knn_ind]=KNN(gnd3,k);
        for i=1:Ndata
            idx=find(gnd3(i,:)>0);
            knn_gnd=gnd3(knn_ind(i,:),:);
            for j=1:length(idx)
                idx_j=idx(j);
                r_score(i,idx_j)= WV(knn_gnd,idx_j,lin_sim);
            end
        end
        [~,ind]=sort(r_score,2,'descend');
        newgnd=gnd3;
        for ii=1:Ndata
            Idx=find(newgnd(ii,:)>0);
            noise=ind(ii,1:num_perprotein_noise(ii));
            childnoise=childGOs(noise);
            newgnd(ii,noise)=0;
            for jj=1:length(childnoise)               % set noisy annotations we identified and their descendants to 0
                childidxp=getGOIdx(childnoise{jj},selGOs);
                childIdx=intersect(childidxp,Idx);
                newgnd(ii,childIdx)=0;            
            end
        end
    
    %% compute precision, recall and f1-measure
    Z =newgnd(1:Ndata,:);
    [tp1,per_pre1,per_re1,per_f11]=PRF(gnd3,Y,Z);
    tp(num_noise,run,:)=tp1;
    per_pre(num_noise,run,:)=per_pre1;
    per_re(num_noise,run,:)=per_re1;
    per_f1(num_noise,run,:)=per_f11;
    Macropre(num_noise,run)=sum(per_pre1)/Ndata;
    Macrore(num_noise,run)=sum(per_re1)/Ndata;
    Macrof1(num_noise,run)=sum(per_f11)/Ndata;
    fprintf('==NoisyGOA1 Run=%d, num_noise=%d, used=%f seconds\n',run,num_noise,etime(clock,t0));
    fprintf('Macropre=%-10.4f,Macrore=%-10.4f, Macrof1=%-10.4f \n\n',...
        Macropre(num_noise,run),Macrore(num_noise,run),Macrof1(num_noise,run));
end
end

Maprecisions=sum(Macropre,2)/round;
Marecall=sum(Macrore,2)/round;
Maf1measure=sum(Macrof1,2)/round;

prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall£¬Maf1';
precision=cell(10,1);
precision{1}=tp;
precision{2}=per_pre;
precision{3}=per_re;
precision{4}=per_f1;
precision{5}=Macropre;
precision{6}=Macrore;
precision{7}=Macrof1;
precision{8}=Maprecisions;
precision{9}=Marecall;
precision{10}=Maf1measure;

stds=cell(10,1);
stds{8}=std(Macropre,0,2);
stds{9}=std(Macrore,0,2);
stds{10}=std(Macrof1,0,2);

if rootGO==5575
    evalstr=['save results',filesep,dataset, '_NoisyGOA1_addNoise1_cc.mat precision stds prec_seq'];
end
if rootGO==3674
    evalstr=['save results',filesep,dataset, '_NoisyGOA1_addNoise1_mf.mat precision stds prec_seq'];
end
if rootGO==8150
    evalstr=['save results',filesep,dataset, '_NoisyGOA1_addNoise1_bp.mat precision stds prec_seq'];
end
eval(evalstr);

fprintf('\n =====finish NoisyGOA1 time=%s\n',datestr(now));

