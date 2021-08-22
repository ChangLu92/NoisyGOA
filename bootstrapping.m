function [maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_miF1score]=bootstrapping(tp,per_pre,per_re,per_f1,num_candidate,num_perprotein_noise)

Ndata=length(per_pre);
ratio=0.85;
time=500;
data=length(find(num_perprotein_noise>0));
% data=229;  %%bp
% data=114;  %%cc
% data=143;  %%mf
selNdata=round(Ndata*ratio);
per_maprecision=zeros(time,selNdata);
per_marecall=zeros(time,selNdata);
per_maF1score=zeros(time,selNdata);
ave_miprecision=zeros(time,1);
ave_mirecall=zeros(time,1);
ave_miF1score=zeros(time,1);

for run=1:time
    sortprotein=randperm(Ndata);
    protein_idx=sortprotein(1:selNdata);
    per_maprecision(run,:)=per_pre(protein_idx);
    per_marecall(run,:)=per_re(protein_idx);  
    per_maF1score(run,:)=per_f1(protein_idx);
    ave_miprecision(run,:)=sum(tp(protein_idx))/sum(num_candidate(protein_idx));
    ave_mirecall(run,:)=sum(tp(protein_idx))/sum(num_perprotein_noise(protein_idx)); 
    ave_miF1score(run,:)=2*ave_miprecision(run,:)*ave_mirecall(run,:)/(ave_miprecision(run,:)+ave_mirecall(run,:));
end
 
ave_maprecision=sum(per_maprecision,2)/data;  
ave_marecall=sum(per_marecall,2)/data;
ave_mafvalue=sum(per_maF1score,2)/data;

maprecisions=sum(ave_maprecision)/time;  
marecalls=sum(ave_marecall)/time;
mafvalue=sum(ave_mafvalue)/time;
miprecisions=sum(ave_miprecision)/time;  
mirecalls=sum(ave_mirecall)/time;
mifvalue=sum(ave_miF1score)/time;

% precision{8}=precisions;
% precision{9}=recalls;
% precision{10}=fvalue;
% precision{11}=ave_precision;
% precision{12}=ave_recall;
% precision{13}=ave_fvalue;
% 
% stds=cell(10,1);
% stds{8}=std(ave_precision,0,1);
% stds{9}=std(ave_recall,0,1);
% stds{10}=std(ave_fvalue,0,1);
