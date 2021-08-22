function [tp,per_pre,per_re,per_f1]=PRF(gnd,Y,Z)
% % compute precision, recall and f1-measure.
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   gnd: the protein-function association matrix after adding noisy annotations 
%   Y: the original protein-function association matrix 
%   Z: the protein-function association matrix after dislodging noisy annotations
%%%%%%%%%%
[Ndata,~]=size(gnd); 
tp=zeros(Ndata,1);
per_pre=zeros(1,Ndata);
per_re=zeros(1,Ndata);
per_f1=zeros(1,Ndata);
for ii=1:size(Y,1)
    noise=find(gnd(ii,:)-Y(ii,:)>0);
    prediction=find((gnd(ii,:)-Z(ii,:))>0);
    if ~isempty(intersect(noise,prediction))
       tp(ii)=length(intersect(noise,prediction));
        per_pre(ii)=tp(ii)/length(prediction);
        per_re(ii)=tp(ii)/length(noise);
        per_f1(ii)=2*per_pre(ii)*per_re(ii)/(per_pre(ii)+per_re(ii));
    end
end