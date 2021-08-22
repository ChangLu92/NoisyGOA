function [] = TestprobNR_hr(dataseth,datasetr,goObj,GOs)
%NR 蛋白质噪声标记剔除
%   dataset是蛋白质功能数据库
%   goObj是GO数据库
%   ccGOs1是GO数据库中cc分支
fprintf('start %s  at %s\n,==Method:%s==',dataseth, datestr(now),'TestprobNR_hr');
load(dataseth);
load(datasetr);

selGOs=GOs;
size_go=length(selGOs);
% lambda=0.5;


if size_go==3765
    gnd1=hGO.ccLabels;
    gnd2=rGO.ccLabels;
    rootGO=5575;   %ccroot
end

if size_go==9878
    gnd1=hGO.mfLabels;
    gnd2=rGO.mfLabels;
    rootGO=3674;   %mfroot
end

if size_go==27342
    gnd1=hGO.bpLabels;
    gnd2=rGO.bpLabels;
    rootGO=8150;  %bproot
end
 
gnd_h=gnd1;
gnd_r=gnd2;

%only test on annotated proteins
index=find(sum(gnd_h,2)==0);
gnd_h(index,:)=[];
gnd_r(index,:)=[];
% [Ndata,Nfun]=size(gnd_h);

minT=1;% the minimum size of member proteins
% maxT=300;
fun_stat_h=sum(gnd_h,1);
fun_stat_r=sum(gnd_r,1);
sel_funh_idx=find(fun_stat_h>=minT);
sel_funr_idx=find(fun_stat_r>=minT);
sel_fun_idx=union(sel_funh_idx,sel_funr_idx);
% sel_fun_idx=intersect(sel_funh_idx,sel_funr_idx);
selGOs=GOs(sel_fun_idx);
gndh=gnd_h(:,sel_fun_idx);
[Ndata, Nfun]=size(gndh);
rootidx=getGOIdx(rootGO,selGOs);

gnd_r=gnd_r(:,sel_fun_idx);
gnd_h=gnd_h(:,sel_fun_idx);

num_perprotein_noise=zeros(Ndata,1); %the number of noisy annotations of each protein
gnd=gnd_r-gnd_h;
sub_goObj=getSelGoObj(selGOs,goObj);%filter the goObj to speedup computation
% DirectchildGOs=getDirectChildGOs(selGOs, selGOs,sub_goObj);
childGOs=getChildGOs(selGOs, selGOs,sub_goObj);
% DirectparGOs=getDirectParentGOs(selGOs, selGOs,sub_goObj);
% Depth  = getSelGOsDepth(selGOs,sub_goObj, rootGO);

%%  count the number of noisy annotations
for ii=1:Ndata
    idx=find(gnd(ii,:)==-1);
    noiseidx{ii}=idx;
    num_perprotein_noise(ii)=length(idx);
end

Y=gnd_r;
gnd3=gnd_h; 

t0=clock;

 %%计算RS正确率，选择所有分数靠前的噪声数量的叶子节点
[newGF, reRatio,preRatio]= probNR(gnd3,selGOs, num_perprotein_noise,noiseidx,childGOs); 
% Fscore=2*preRatio*reRatio/(preRatio+reRatio);
Z =newGF(1:Ndata,:);

[tp,per_pre,per_re,per_f1,Miprecisions,Mirecall,num_candidate]=PRF1(gnd_h,Y,Z);
data=length(find(num_perprotein_noise>0));
[maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_mifvalue]=bootstrapping(tp,per_pre,per_re,per_f1,num_candidate,num_perprotein_noise);
prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall,Maf1,Miprecisions,Mirecall，Mif1,maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_mifvalue,num_perprotein_noise,num_candidate';

precision=cell(30,1);
precision{1}=tp;
precision{2}=per_pre;
precision{3}=per_re;
precision{4}=per_f1;
precision{5}=sum(per_pre)/data;
precision{6}=sum(per_re)/data;
precision{7}=sum(per_f1)/data;
precision{8}=Miprecisions;
precision{9}=Mirecall;
precision{10}=2*Miprecisions*Mirecall/(Miprecisions+Mirecall);
precision{11}=maprecisions;
precision{12}=marecalls;
precision{13}=mafvalue;
precision{14}=miprecisions;
precision{15}=mirecalls;
precision{16}=mifvalue;
precision{17}=ave_maprecision;
precision{18}=ave_marecall;
precision{19}=ave_mafvalue;
precision{20}=ave_miprecision;
precision{21}=ave_mirecall;
precision{22}=ave_mifvalue;
precision{23}=num_perprotein_noise;
precision{24}=num_candidate;


% 
stds=cell(30,1);
stds{11}=std(ave_maprecision,0,1);
stds{12}=std(ave_marecall,0,1);
stds{13}=std(ave_mafvalue,0,1);
stds{14}=std(ave_miprecision,0,1);
stds{15}=std(ave_mirecall,0,1);
stds{16}=std(ave_mifvalue,0,1);


if rootGO==5575
evalstr=['save results',filesep,dataseth, '_probNR_cc.mat precision stds prec_seq'];
end
if rootGO==3674
evalstr=['save results',filesep,dataseth, '_probNR_mf.mat precision stds prec_seq'];
end
if rootGO==8150
evalstr=['save results',filesep,dataseth, '_probNR_bp.mat precision stds prec_seq'];
end
eval(evalstr);

     
fprintf('\n =====finish probNR time=%s\n',datestr(now))
%end

