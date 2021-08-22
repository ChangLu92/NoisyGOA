function [lin_sim] = getLinStructSim(selGOs,childGOs,parGOs, nTerms) 
%% compute the Lin's structure similarity
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   selGOs: the selected GOs vector
%   childGOs: the children GOs set of the selected GOs
%   parGOs: the ancestors GOs set of the selected GOs
%   nTerms:  the number of GO terms involved
%%%%%%%%%%
[Nfun] = length(selGOs);
structICs = compStructIC(childGOs,  selGOs, nTerms);    %calculate the struct-IC of the selected GO terms  
%% calculate the lin_sim via the structed IC       
lin_sim=zeros(Nfun,Nfun);
goids=num2goid(selGOs);
for ii = 1:length(selGOs)
    for jj = ii+1:length(selGOs)
        ancestors=intersect(parGOs{ii}, parGOs{jj});
        p_idx=getGOIdx(ancestors,selGOs);
        p_idx(p_idx==0)=[];
        if ~isempty(p_idx) 
            [pca,min_idx]=max(structICs(p_idx));%the most specific common ancestor
            sim=2*pca/(structICs(ii)+structICs(jj));
            lin_sim(ii,jj)=sim;
            lin_sim(jj,ii)=sim;
        else
            fprintf('no shared ancestor for %s and %s\n', goids{ii}, goids{jj});
        end
    end
    lin_sim(ii,ii)=1;
end


