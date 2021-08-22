
load('GeneOntointer2016102920170331.mat');

size_go=length(ccGOs);
childGOs=getChildGOs(ccGOs, ccGOs,goObj);
parGOs=getParentGOs(ccGOs, ccGOs,goObj);
lin_sim_cc= getLinStructSim(ccGOs,childGOs,parGOs, size_go);

size_go=length(mfGOs);
childGOs = getChildGOs(mfGOs, mfGOs,goObj);
parGOs = getParentGOs(mfGOs, mfGOs,goObj);
lin_sim_mf= getLinStructSim(mfGOs,childGOs,parGOs, size_go);

size_go=length(bpGOs);
childGOs=getChildGOs(bpGOs, bpGOs,goObj);
parGOs=getParentGOs(bpGOs, bpGOs,goObj);
lin_sim_bp= getLinStructSim(bpGOs,childGOs,parGOs, size_go);

save lin_sim_struct_cc.mat lin_sim_cc;
save lin_sim_struct_mf.mat lin_sim_mf;
save lin_sim_struct_bp.mat lin_sim_bp;
