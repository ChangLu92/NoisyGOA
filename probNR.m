function [ newGF ,recall,precision ] = probNR(gnd, selGOs, num_perprotein_noise,noise,childGOs)
%按照功能概率选择节点删除，不是叶子节点的，将孩子节点全部删除

allrinoise=0;
num_candidate=0;
[Ndata,Nfun]=size(gnd);
probfun=zeros(1,Nfun);
fun_count=sum(gnd); %每一个term被多少蛋白质标注
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
        childnoise=childGOs(noiseIdx);  %找p个噪声节点的孩子节点
        childIdx=[];
        for jj=1:length(childnoise)
            childidxj=getGOIdx(childnoise{jj},selGOs); %孩子节点的序号
            childidx=intersect(childidxj,ind_i);    %蛋白质ii中存在的孩子节点
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
     recall=allrinoise/sum(num_perprotein_noise);   %recall:正确的噪声/全部实际噪声
     precision=allrinoise/num_candidate;            %precision：正确的噪声/全部找到的噪声
fprintf('probNR recall= %d , precision= %d \n',recall,precision);
fprintf('probNR is over! \n');
end