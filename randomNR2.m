function [ gnd2 , recall, precision] = randomNR2(gnd2, Ndata,childGOs, selGOs, num_perprotein_noise,noise)
%%随机寻找num_perprotein_noise个节点,不是叶子节点将其所有孩子节点设置为噪声
% 
allrinoise=0;
num_candidate=0;
for ii=1:Ndata
    Idx=find(gnd2(ii,:)>0);
    p=randperm(length(Idx));
    p=p(1:num_perprotein_noise(ii)); 
    %ran_noise=selGOs(Idx(p));    %随机找出p个噪声
    %rannoiseIdx=getGOIdx(ran_noise,selGOs); 
    rannoiseIdx=Idx(p);
    childnoise=childGOs(rannoiseIdx);  %找p个噪声节点的孩子节点
    childIdx=[];
    for jj=1:length(childnoise)
        childidxj=getGOIdx(childnoise{jj},selGOs); %孩子节点的序号
        childidx=intersect(childidxj,Idx);    %蛋白质ii中存在的孩子节点
        if ~isempty(childidx)
           childIdx=[childIdx;childidx];
        end 
    end
    childIdx=childIdx';
%     [a,b]=size(rannoiseIdx);
%     [c,d]=size(childIdx);
    %fprintf('a= %d ,b= %d ,c= %d ,d= %d \n',a,b,c,d);
    candidate=[rannoiseIdx,childIdx];  %%%%%这里原本的逗号改成了分号
    candidate=unique(candidate);
    rightnoise=intersect(candidate,noise{ii});
    gnd2(ii,candidate)=0;
    num_candidate=length(candidate)+num_candidate;
    allrinoise=length(rightnoise)+allrinoise;
end
     recall=allrinoise/sum(num_perprotein_noise);   %recall:正确的噪声/全部实际噪声
     precision=allrinoise/num_candidate;            %precision：正确的噪声/全部找到的噪声
     fprintf('rNR recall= %d , precision= %d',recall,precision);
     fprintf('randomNR2 is over!');
end
    
    