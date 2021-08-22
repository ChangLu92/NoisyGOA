function [ gnd2 , recall, precision] = randomNR2(gnd2, Ndata,childGOs, selGOs, num_perprotein_noise,noise)
%%���Ѱ��num_perprotein_noise���ڵ�,����Ҷ�ӽڵ㽫�����к��ӽڵ�����Ϊ����
% 
allrinoise=0;
num_candidate=0;
for ii=1:Ndata
    Idx=find(gnd2(ii,:)>0);
    p=randperm(length(Idx));
    p=p(1:num_perprotein_noise(ii)); 
    %ran_noise=selGOs(Idx(p));    %����ҳ�p������
    %rannoiseIdx=getGOIdx(ran_noise,selGOs); 
    rannoiseIdx=Idx(p);
    childnoise=childGOs(rannoiseIdx);  %��p�������ڵ�ĺ��ӽڵ�
    childIdx=[];
    for jj=1:length(childnoise)
        childidxj=getGOIdx(childnoise{jj},selGOs); %���ӽڵ�����
        childidx=intersect(childidxj,Idx);    %������ii�д��ڵĺ��ӽڵ�
        if ~isempty(childidx)
           childIdx=[childIdx;childidx];
        end 
    end
    childIdx=childIdx';
%     [a,b]=size(rannoiseIdx);
%     [c,d]=size(childIdx);
    %fprintf('a= %d ,b= %d ,c= %d ,d= %d \n',a,b,c,d);
    candidate=[rannoiseIdx,childIdx];  %%%%%����ԭ���Ķ��Ÿĳ��˷ֺ�
    candidate=unique(candidate);
    rightnoise=intersect(candidate,noise{ii});
    gnd2(ii,candidate)=0;
    num_candidate=length(candidate)+num_candidate;
    allrinoise=length(rightnoise)+allrinoise;
end
     recall=allrinoise/sum(num_perprotein_noise);   %recall:��ȷ������/ȫ��ʵ������
     precision=allrinoise/num_candidate;            %precision����ȷ������/ȫ���ҵ�������
     fprintf('rNR recall= %d , precision= %d',recall,precision);
     fprintf('randomNR2 is over!');
end
    
    