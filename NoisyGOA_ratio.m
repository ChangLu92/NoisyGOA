function [] = NoisyGOA_ratio(dataset,goObj,GOs)
% Identifying noisy gene ontology annotations using taxonomic and semantic
% similarity by adding fixed retio of annotations
% Chang Lu, Jun Wang, Guoxian Yu, College of Computer and Information
% Science, Southwest University. Contact gxyu@swu.edu.cn, lucifer1992@email.swu.edu.cn
%%%%%%%%%%
%   goObj: gene ontology read from GO file (obo format)
%   GOs:  GO term id for each column of the protein-function association matrix
%   dataset: Name of the dataset, 'YeastGOA' means Yeast proteins annotated with terms.
%%%%%%%%%%
fprintf('start %s  at %s\n,==Method:%s==',dataset, datestr(now),'NoisyGOA2');
load(dataset);
selGOs=GOs;
size_go=length(selGOs);

if size_go==3865
    gnd=GO.ccLabels;
    rootGO=5575;   %ccroot
    load('lin_sim_struct_cc.mat', 'lin_sim');
end

if size_go==9982
    gnd=GO.mfLabels;
    rootGO=3674;   %mfroot
    load('funsim_lin_struct_mf.mat', 'lin_sim');
end

if size_go==28195
    gnd=GO.bpLabels;
    rootGO=8150;  %bproot
    load('funsim_lin_struct_bp.mat', 'lin_sim');
end

percentage=0.05:0.05:0.3; % Ratio of noisy annotations
k=5;              % number of nearest neighbourhood proteins
round=5;          % time of the independent repeated experiment

%only test on annotated proteins
gnd2=full(gnd);
index=find(sum(gnd2,2)==0);
gnd2(index,:)=[];
[Ndata,Nfun]=size(gnd2);

tp=zeros(length(percentage),round,Ndata); %the number of correctly predicted noisy annotations of every protein
per_pre=zeros(length(percentage),round,Ndata); % precision of identifying noisy annotations of every protein
per_re=zeros(length(percentage),round,Ndata);  % recall of identifying noisy annotations of every protein
per_f1=zeros(length(percentage),round,Ndata);  % f1-measure of identifying noisy annotations of every protein
Macropre=zeros(length(percentage),round);    % average precision
Macrore=zeros(length(percentage),round);     % average recall
Macrof1=zeros(length(percentage),round);     % average f1-measure

childGOs=getChildGOs(selGOs, selGOs,goObj);
DirectchildGOs=getDirectChildGOs(selGOs, selGOs,goObj);    %寻找直系孩子节点

for percentage=0.05:0.05:0.3
    num=percentage*20;
    for run=1:round
        %% add noisy annotations
         num_noise_Mat=round(sum(gnd2,2)*percentage);
        [gnd3,num_perprotein_noise] = addNoise2(gnd2,selGOs,num_noise_Mat,DirectchildGOs);
        fprintf('==finish addNoise!\n');
        
        %% Identifying Noisy Gene Ontology Annotations
        Y=gnd2;
        gnd3=full(gnd3);
        t0=clock;
        r_score=zeros(Ndata,Nfun);
        [~,p_knind]=KNN(gnd3,k);
        for i=1:Ndata
            idx=find(gnd3(i,:)>0);
            knn_gnd=gnd3(p_knind(i,:),:);
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
            childIdx=[];
            for jj=1:length(childnoise)
                childidxp=getGOIdx(childnoise{jj},selGOs);
                childidx=intersect(childidxp,Idx);
                if ~isempty(childidx)
                    childIdx=[childIdx;childidx];
                end
            end
            childIdx=childIdx';
            candidate=[noise,childIdx];
            candidate=unique(candidate);
            newgnd(ii,candidate)=0;
        end
        
        %% compute precision, recall and f1-measure
        Z =newgnd(1:Ndata,:);  
        [tp1,per_pre1,per_re1,per_f11]=PRF(gnd3,Y,Z);
        tp(num,run,:)=tp1;
        per_pre(num,run,:)=per_pre1;
        per_re(num,run,:)=per_re1;
        per_f1(num,run,:)=per_f11;
        Macropre(num,run)=sum(per_pre1)/Ndata;
        Macrore(num,run)=sum(per_re1)/Ndata;
        Macrof1(num,run)=sum(per_f11)/Ndata;
        fprintf('==NoisyGOA2 Run=%d, percentage=%-10.4f, used=%f seconds\n',run,percentage,etime(clock,t0));
        fprintf('Macropre=%-10.4f,Macrore=%-10.4f, Macrof1=%-10.4f \n\n',...
            Macropre(num,run),Macrore(num,run),Macrof1(num,run));
    end
end

Maprecisions=sum(Macropre,2)/round;
Marecall=sum(Macrore,2)/round;
Maf1measure=sum(Macrof1,2)/round;

prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall，Maf1';
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
    evalstr=['save results',filesep,dataset, '_NoisyGOA2_addNoise2_cc.mat precision stds prec_seq'];
end
if rootGO==3674
    evalstr=['save results',filesep,dataset, '_NoisyGOA2_addNoise2_mf.mat precision stds prec_seq'];
end
if rootGO==8150
    evalstr=['save results',filesep,dataset, '_NoisyGOA2_addNoise2_bp.mat precision stds prec_seq'];
end
eval(evalstr);

fprintf('\n =====finish NoisyGOA2 time=%s\n',datestr(now));

