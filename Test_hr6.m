% clear all;
datapath=[pwd,filesep,'data',filesep];%pwd is the current work directory
addpath(datapath);                    %pwd 当前路径；filesep 当前平台（win or linux）‘/’or '\'
% datapath=[pwd,filesep,'funs1',filesep];%pwd is the current work directory
% addpath(datapath); 
datapath=[pwd,filesep,'funs',filesep];%pwd is the current work directory
addpath(datapath);

load('GeneOntointer201505250912.mat');
% EC={'ZebrafishGOAinternotH20150105ZebrafishGOAinternotR20150526ECs.mat','MouseGOAinternotH20150105MouseGOAinternotR20150602ECs.mat','ArabidopsisGOAinternotH20150105ArabidopsisGOAinternotR20150330ECs_sum'};
 datasets={'ZebrafishGOAinternotH20150526';'ZebrafishGOAinternotR20150914'};
  datasets2={'MouseGOAinternotH20150602';'MouseGOAinternotR20150914'};
% datasets3={'ArabidopsisGOAinternotH20150526';'ArabidopsisGOAinternotR20150914'};
% times=zeros(length(datasets),5);
goObj_h=goObj20150525;
GOs{1}=ccGOs;
GOs{2}=mfGOs;
GOs{3}=bpGOs;

for ii=1:3
%  a=0.2;
%  w=0.5;
%  lambda=0.5;
%  for  a=0:0.1:0.9
%   NoisyGOA_NtNEC( datasets{1},datasets{2},EC{1},goObj_h,GOs{ii},a);
% NoisyGOA_NtNEC( datasets2{1},datasets2{2},EC{2},goObj_h,GOs{ii},a);
% NoisyGOA_NtNEC( datasets3{1},datasets3{2},EC{3},goObj_h,GOs{ii},a);
%       NoisyGOA_hrEC1( datasets{1},datasets{2},EC{1},goObj_h,GOs{ii},a);
%      NoisyGOA_hrEC1( datasets2{1},datasets2{2},EC{2},goObj_h,GOs{ii},a);
%      NoisyGOA_hrEC1( datasets3{1},datasets3{2},EC{3},goObj_h,GOs{ii},a);
% TestranNR2_hr(datasets{1},datasets{2},goObj_h,GOs{ii});
% TestranNR2_hr(datasets2{1},datasets2{2},goObj_h,GOs{ii});
TestprobNR_hr(datasets{1},datasets{2},goObj_h,GOs{ii});
TestprobNR_hr(datasets2{1},datasets2{2},goObj_h,GOs{ii});
%         NoisyGOA_ECSP5norm(datasets{1},datasets{2},EC{1}, goObj_h,GOs{ii},a,w,lambda);
%     NoisyGOA_ECSP5norm(datasets2{1},datasets2{2},EC{2}, goObj_h,GOs{ii},a,w,lambda);
%      NoisyGOA_ECSP5norm2(datasets3{1},datasets3{2},EC{3}, goObj_h,GOs{ii},a,w,lambda);
%   NoisyGOA_SP1(datasets{1},datasets{2},goObj_h,GOs{ii},lambda);
%   NoisyGOA_SP1(datasets2{1},datasets2{2},goObj_h,GOs{ii},lambda);
%   NoisyGOA_SP1(datasets3{1},datasets3{2},goObj_h,GOs{ii},lambda);
%    NoisyGOA_ECSP5norm(datasets{1},datasets{2},EC{1}, goObj_h,GOs{ii},a,w);
%     NoisyGOA_ECSP5norm(datasets2{1},datasets2{2},EC{2}, goObj_h,GOs{ii},a,w);
%      NoisyGOA_ECSP5norm(datasets3{1},datasets3{2},EC{3}, goObj_h,GOs{ii},a,w);
%    NoisyGOA_EC2015(datasets{1},datasets{2},EC{1},goObj_h,GOs{ii});
%    NoisyGOA_EC2015(datasets2{1},datasets2{2},EC{2},goObj_h,GOs{ii});
%    NoisyGOA_EC2015(datasets3{1},datasets3{2},EC{3},goObj_h,GOs{ii});
   %times(ii,1)=etime(clock,t0);
%  end
end


     