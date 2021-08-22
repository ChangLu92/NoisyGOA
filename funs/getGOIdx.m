function [ goIdx ] = getGOIdx(go, GOs)
%  get the GO's the index (according to the integer label) in the selGOs,
%  all input in the numeric
%  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
%   Southwest University.
%   version 1.0 date:2014-02-26
%%%%%%%%%%
%   go: GO term id
%   GOs: the selected GOs vector
%%%%%%%%%%
goIdx=zeros(length(go),1);
for ii=1:length(go)
    goii=go(ii);
    [idx]=find(GOs==goii,1); 
    if(~isempty(idx))
        goIdx(ii)=idx;
    end
end
