function [gnd]=NAI(gnd,num_perprotein_noise,selGOs,childGOs,ind)
% % Noisy Annotations Identification. select designated terms and their descendants as noisy annotations.
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
for ii=1:Ndata
    Idx=find(gnd(ii,:)>0);
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
    gnd(ii,candidate)=0;
end
fprintf('NAI is over!');
end