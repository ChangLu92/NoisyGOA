function [ parGOs ] = getParentGOs(childGO, selGOs, goObj)
% % %get the parentGOs of childGO(without prefix 'GO:', all input in the numeric form).
% % %  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2014-02-26
% childGO=[5976;8150];
% selGOs=[5976;8150;8152;71704];
% yeastStruct=goannotread('YeastGO.txt','Aspect','P','Fields',{'GOid'});
%goObj=geneont('File','go20140201.obo');

num_child=size(childGO,1);
parGOs=cell(size(childGO,1),1);
for ii=1:num_child
    child=childGO(ii);
    parGO=getancestors(goObj,child,'Exclude',true);
    parGOs{ii}=intersect(parGO,selGOs);
end

%parGoIDs=cellstr(num2goid(parGoNums));
%parNodeNames=struct2cell();

% parNodeNames=[];
% pattern='.';
% childNodeName=char(childNodeName);
% if(length(childNodeName)>2)% 2nd parent nodes
%     k=strfind(childNodeName,pattern);
%     levels=size(k,2);
%     %levels=length(k_ch);
%     parNodeNames=cell(levels,1);
%     jj=1;
%     for ii=levels:-1:1%from the parents to the grandparents
%         parent=childNodeName(1:k(ii)-1);
%         k2=intersect(selGONames,parent);
%         sel=0;
%         k2_length=size(k2,1);
%        % k2_length=length(k2_ch);
%         if(k2_length>0)
%             sel=1;
%         end
%         if sel
%             parNodeNames(jj)=cellstr(parent);
%             jj=jj+1;
%         end
%     end
%     parNodeNames(jj:levels)=[];
%end