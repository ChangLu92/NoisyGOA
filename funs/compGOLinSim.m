 function [ sim ] = compGOLinSim(go1, go2, parGO1, parGO2,selGOs, selGOProbs)
% compute the Lin's similarity. ICML98, Dekang Lin's paper
%  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
%   Southwest University.
%   version 1.0 date:2014-02-26
 %goObj=geneont('File','go20140201.obo');
 %goObj=geneont('File','GO20140201.obo');
% selGOProbs=[0.1;0.2;0.3;0.4;0.5];
% selGOs=[5975;5976;8150;8152;71704];
% go1=5975;
% go2=5976;

 %parGO1=getParentGOs(go1,selGOs,goObj);
% parGO2=getParentGOs(go2,selGOs,goObj);

parGOs=intersect(parGO1,parGO2);
sim=0;
if ~isempty(parGOs)%have parents
    par_idx=getGOIdx(parGOs,selGOs);
    if ~isempty(par_idx)
        [pca,min_idx]=min(selGOProbs(par_idx));%the most specific common ancestor
        go_idx1=getGOIdx(go1,selGOs);
        go_idx2=getGOIdx(go2,selGOs);
        p1=selGOProbs(go_idx1);
        p2=selGOProbs(go_idx2);
        sim=2*log(pca)/(log(p1)+log(p2));
    end
else
    sim=0;
end

if isnan(sim)
    sim=0;
end
if isinf(sim)
	sim=0;
end

