function [  r_score ] = WV( knn_y,y,sim_funs)
% % compute a score leveraging the semantic similarity and the taxonomic similarity between two terms.
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   knn_y: k nearest neighbourhood protein of protein y
%   y: a annotation of protein y
%   sim_funs: taxonomic similarity between two terms 
%%%%%%%%%%
[k,n]=size(knn_y); 
sim=repmat(sim_funs(y,:),k,1);
sim_y=knn_y.*sim;
r_score=sum(max(sim_y')); 
if r_score==0
   r_score=eps;
end
r_score=1/r_score; 
end

