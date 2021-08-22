 function [ childGOs ] = getDirectChildGOs(parGO, selGOs, goObj)
% % %  get the parentGOs of direct childGO.
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   parGO:  the parents GOs set of the selected GOs
%   selGOs: the selected GOs vector
%   goObj:  gene ontology read from GO file (obo format)
%%%%%%%%%%
num_child=size(parGO,1);
childGOs=cell(num_child,1);
for ii=1:num_child
    go=parGO(ii);
    childGO=getdescendants(goObj, go, 'Depth', 1,'Exclude',true);
    childGOs{ii}=intersect(childGO,selGOs);
end

