datapath=[pwd,filesep,'data',filesep];%pwd is the current work directory
addpath(datapath);                    
datapath=[pwd,filesep,'funs',filesep];
addpath(datapath); 

load('GeneOnto.mat');% load GO file 
%dataset='YeastGOA'; 
dataseth ='YeastGOA2h';%historical GOA, archieved date 2015-05-26
datasetr='YeastGOA2r';%recent GOA, archieved date 2015-12-07
funGOs=ccGOs; %GO term numeric id 

%  NoisyGOA_num(dataseth,goObj,funGOs);
%  NoisyGOA_ratio(dataseth,goObj,funGOs);
%  NoisyGOA_hr(dataseth,datasetr,goObj,funGOs);
% LeastVote(dataseth,goObj,funGOs);
NoisyGOAlin_hr(dataseth,datasetr,goObj,funGOs);

