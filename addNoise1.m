function [ gnd ,num_perprotein_noise] = addNoise_num(gnd,selGOs,num_noise,DirectchildselGOs)
% % adding fixed number of annotations of each protein.
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   DirectchildselGOs: the direct children GOs set of the selected GOs
%   selGOs: the selected GOs vector
%   num_noise: the number of noisy annotations we prepare to add
%   num_perprotein_noise: the number of noisy annotations we added actually
%   gnd: a protein-function association matrix
%%%%%%%%%%
[Ndata,~]=size(gnd); 
num_perprotein_noise=zeros(Ndata,1); %number of noise annotations actually we added
for ii=1:Ndata
    num_n=num_noise;    
    if ~isempty(find(gnd(ii,:)>0, 1));
        lengthnoise=0;
        while num_n>0            
            idx=find(gnd(ii,:)>0);
            gos=selGOs(idx);              
            child=DirectchildselGOs(idx);    
            diffGOs=[];
            for jj=1:length(idx)         %find all the direct children GOs of this protein
                childGO=child{jj};       
                differentGOs=setdiff(childGO,gos);
                diffGOs=[diffGOs;differentGOs];
            end
            diffGOs=unique(diffGOs);
            ran=randperm(length(diffGOs));
            if (ran>0)
                ran=ran(1);
                noiseGOs=diffGOs(ran);
                noise1=getGOIdx(noiseGOs,selGOs);
                gnd(ii,noise1)=1;
                lengthnoise=lengthnoise+1;
            end
            num_n=num_n-1;
        end
        num_perprotein_noise(ii)=num_perprotein_noise(ii)+lengthnoise;
    end
end
