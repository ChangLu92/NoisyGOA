function [ newGF ,recall,precision ] = probNR(gnd, selGOs, num_perprotein_noise,noise,childGOs)
%���չ��ܸ���ѡ��ڵ�ɾ��������Ҷ�ӽڵ�ģ������ӽڵ�ȫ��ɾ��

allrinoise=0;
num_candidate=0;
[Ndata,Nfun]=size(gnd);
probfun=zeros(1,Nfun);
fun_count=sum(gnd); %ÿһ��term�����ٵ����ʱ�ע
fun_count_idx=find(fun_count~=0);
probfun(fun_count_idx)=1./fun_count(fun_count_idx);
freq= repmat(probfun,Ndata,1);
newgnd=gnd.*freq;
[value,ind]=sort(newgnd,2,'descend'); 

for ii=1:Ndata
    if ~isempty(num_perprotein_noise(ii))
        idx=find(value(ii,:)~=0);
        ind_i=ind(ii,idx);
        noiseIdx=ind_i(1:num_perprotein_noise(ii));
        noiseIdx=noiseIdx';
        childnoise=childGOs(noiseIdx);  %��p�������ڵ�ĺ��ӽڵ�
        childIdx=[];
        for jj=1:length(childnoise)
            childidxj=getGOIdx(childnoise{jj},selGOs); %���ӽڵ�����
            childidx=intersect(childidxj,ind_i);    %������ii�д��ڵĺ��ӽڵ�
        if ~isempty(childidx)
           childIdx=[childIdx;childidx];
%            childIdx=[childIdx,childidx];
        end 
        end
   % childIdx=childIdx';
    candidate=[noiseIdx;childIdx];
    candidate=unique(candidate);
    end
    rightnoise=intersect(candidate,noise{ii});
    gnd(ii,candidate)=0;
    num_candidate=length(candidate)+num_candidate;
    allrinoise=length(rightnoise)+allrinoise;
end
newGF=gnd;
     recall=allrinoise/sum(num_perprotein_noise);   %recall:��ȷ������/ȫ��ʵ������
     precision=allrinoise/num_candidate;            %precision����ȷ������/ȫ���ҵ�������
fprintf('probNR recall= %d , precision= %d \n',recall,precision);
fprintf('probNR is over! \n');
end