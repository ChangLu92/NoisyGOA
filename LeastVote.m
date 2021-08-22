function [] = LeastVote(dataset,goObj,newccGOs)
%LeastVote 蛋白质噪声标记剔除
%   dataset是蛋白质功能数据库
%   goObj是GO数据库
%   ccGOs1是GO数据库中cc分支
fprintf('start %s  at %s\n,==Method:%s==',dataset, datestr(now),'TestrsvoteNR4');
load(dataset);
ccGOs=newccGOs;
size_go=length(ccGOs);
if size_go==3865
    gnd=hGO.ccLabels;
    rootGO=5575;   %ccroot
end

if size_go==9982
    gnd=GO.mfLabels;
    rootGO=3674;   %mfroot
end

if size_go==28195
    gnd=hGO.bpLabels;
    rootGO=8150;  %bproot
end 

rootIdx=getGOIdx(rootGO,ccGOs);  %在ccGO中找root节点的idx（是第几个）
num_noise=1:5;   %噪声的数量
round=5;
gnd2=full(gnd);
index=find(sum(gnd2,2)==0);
%only test on annotated proteins
gnd2(index,:)=[];
[Ndata,Nfun]=size(gnd2);
% [Ndata,Nfun]=size(gnd);
rsprecision1=zeros(length(num_noise),round);
recall1=zeros(length(num_noise),round);
F1score1=zeros(length(num_noise),round);
tp=zeros(length(num_noise),round,Ndata);
per_pre=zeros(length(num_noise),round,Ndata);
per_re=zeros(length(num_noise),round,Ndata);
per_f1=zeros(length(num_noise),round,Ndata);
Macropre=zeros(length(num_noise),round);
Macrore=zeros(length(num_noise),round);
Macrof1=zeros(length(num_noise),round);

selGOs=ccGOs;

childGOs=getChildGOs(ccGOs, ccGOs,goObj);
DirectchildGOs=getDirectChildGOs(ccGOs, ccGOs,goObj);    %寻找直系孩子节点


for num_noise=3:5
    for run=1:round
        %添加噪声
       [gnd3, num_perprotein_noise]=addNoise1(gnd2,selGOs,num_noise,DirectchildGOs);       
       fprintf('==finish New_addNoise!\n');
        Y=gnd2;
        t0=clock;
        r_score=zeros(Ndata,Nfun);
        [p_knsor,p_knind]=KNN(gnd3,5); %取k为8，p_knsor是排序的值,p_knind是排序的顺序
        for i=1:Ndata
             Idx=find(gnd3(i,:)==0);
            knn_gnd=gnd3(p_knind(i,:),:);
            sum_knn_gnd=sum(knn_gnd);
             sum_knn_gnd(Idx)=0;
            sum_knn_gnd=sum_knn_gnd+gnd3(i,:);
            b=sum_knn_gnd~=0;
            c=1./sum_knn_gnd(b);
            r_score(i,b)=c;
        end
        [~,ind]=sort(r_score,2,'descend');
        newgnd=gnd3;
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

        Z =newgnd(1:Ndata,:);
        [tp1,per_pre1,per_re1,per_f11]=PRF(gnd3,Y,Z);
        tp(num_noise,run,:)=tp1;
        per_pre(num_noise,run,:)=per_pre1;
        per_re(num_noise,run,:)=per_re1;
        per_f1(num_noise,run,:)=per_f11;
        Macropre(num_noise,run)=sum(per_pre1)/Ndata;
        Macrore(num_noise,run)=sum(per_re1)/Ndata;
        Macrof1(num_noise,run)=sum(per_f11)/Ndata;
        fprintf('==deepestNR Run=%d, num_noise=%-10.4f, used=%f seconds\n',run,num_noise,etime(clock,t0));
        fprintf('Macropre=%-10.4f,Macrore=%-10.4f, Macrof1=%-10.4f \n\n',...
        Macropre(num_noise,run),Macrore(num_noise,run),Macrof1(num_noise,run));
    end
end
precisions1=sum(rsprecision1,2)/round;
recalls1=sum(recall1,2)/round;
fvalue1=sum(F1score1,2)/round;
Maprecisions=sum(Macropre,2)/round;
Marecall=sum(Macrore,2)/round;
Maf1=sum(Macrof1,2)/round;
prec_seq='precisions,recall,F1score';
precision=cell(15,1);
precision{1}=precisions1;
precision{2}=recalls1;
precision{3}=fvalue1;
precision{4}=tp;
precision{5}=per_pre;
precision{6}=per_re;
precision{7}=per_f1;
precision{8}=Macropre;
precision{9}=Macrore;
precision{10}=Macrof1;
precision{11}=Maprecisions;
precision{12}=Marecall;
precision{13}=Maf1;

stds=cell(13,1);
stds{1}=std(rsprecision1,0,2);
stds{2}=std(recall1,0,2);
stds{3}=std(F1score1,0,2);
stds{8}=std(Macropre,0,2);
stds{9}=std(Macrore,0,2);
stds{10}=std(Macrof1,0,2);

if rootGO==5575
evalstr=['save results',filesep,dataset, '_rsvoteallNR_newaddno_cck5_prf.mat precision stds prec_seq'];
end
if rootGO==3674
evalstr=['save results',filesep,dataset, '_rsvoteallNR_newaddno_mfk5_prf.mat precision stds prec_seq'];
end
if rootGO==8150
evalstr=['save results',filesep,dataset, '_rsvoteallNR_newaddno_bpk5_prf.mat precision stds prec_seq'];
end
eval(evalstr);


fprintf('\n =====finish LV time=%s\n',datestr(now))
%end

