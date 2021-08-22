function [ knsor,knind ] = KNN(x,k)
%compute k nearest neighbour proteins of x
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   x: a protein vector
%   k: the number of nearest neighbour proteins
%%%%%%%%%%
distance= pdist2(x,x,'cosine'); 
[sorted,index]=sort(distance,2);
knsor=sorted(:,2:k+1);
knind=index(:,2:k+1);
fprintf('KNN is over !\n');


