function [ structICs ] = compStructIC(descendantGOs,  selGOs)
% % compute the structure information content of each selected terms.
% % %  Coded by Guangyuan Fu (fugy@swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 2.0 date:2016-01-26
%%%%%%%%%%
%   descendantGOs: the children GOs set of the selected GOs
%   selGOs: the selected GOs vector
%   nTerms:  the number of GO terms involved
%%%%%%%%%% 
nTerms=length(selGOs);
 structICs=ones(length(selGOs),1);
 for ii=1:length(selGOs)
     if ~isempty(descendantGOs{ii})
        structICs(ii)=(1-log(length(descendantGOs{ii})+1)/log(nTerms));
     else
      structICs(ii)=1;
     end
 end
 structICs( isnan(structICs))=0;
 structICs( isinf(structICs))=0;