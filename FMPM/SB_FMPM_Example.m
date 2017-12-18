clear all;
close all;
fclose all;
addpath('F:\SSPM-V2.1\FMPM\');
workDir='\RealData_NDARTBSS';
PREfile='\RealData_NDAR';
APENfile='.mat';
DataFile=strcat('F:\SSPM-V2.1\FMPM',PREfile,APENfile);
load(DataFile); 
Xdesign=[Xdesign(:,1:2), Xdesign(:,3), Xdesign(:,3).^2];
Zdesign= Xdesign(:,[1,3]);
[NoSetup, Xdesign, Zdesign, Ydesign, scalediffusion]=FMEM_ReadNormalize(Xdesign,Zdesign,Ydesign,Indicator, arclength);
n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
sumri=NoSetup(4);
pz=NoSetup(5);
ExpVar=0.99;

GG=5000;
rng('default');
[mh1,GCVs,vh] = FMPM_LL_wob(NoSetup,Indicator,arclength,Xdesign,Ydesign);
[efitBetas,efitBetas1,InvSigmats,efitYdesign,ResYdesign] = FMPM_LL_wb(NoSetup,arclength,Xdesign,Ydesign,mh1);

[SigmaZeta,SigmaZetad,SigmaXi,SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,ResYdesign);
[mh0] = BiSmooth_LL_bs(NoSetup,arclength,SigmaXi,pz,0);
[SigmaXiP,SigmaXiE] = BiSmooth_LL(NoSetup,arclength,SigmaXi,SigmaXid,pz,mh0,0);
[mh00] = BiSmooth_LL_bs(NoSetup,arclength,SigmaZeta,1,0);
[SigmaZetaP,SigmaZetaE] = BiSmooth_LL(NoSetup,arclength,SigmaZeta, SigmaZetad,1,mh00,0);
SigmaXiE=kron(eye(L0),SigmaXiE);
SigmaZetaE=SigmaZetaE*eye(L0);
[~,~,~,~,~,~,SigmaXiPCA, SigmaZetaPCA] = FMPM_FPCA_Cov(NoSetup,SigmaXiP,SigmaZetaP,SigmaXiE,SigmaZetaE,ExpVar);


[mh4,GCVs,vh,vYdesign, tempv] = FMPM_LL_Ref_wob_SB_T(NoSetup,Indicator, arclength,Xdesign,Zdesign,Ydesign, SigmaXiPCA,SigmaZetaPCA);
[efitBetas,efitBetas1,efitYdesign, ResYdesign,InvSigmats,Tempmat]= FMPM_LL_Ref_wb_SB_T(NoSetup,Indicator, arclength,Xdesign,Zdesign, Ydesign,vYdesign,tempv,SigmaXiPCA,SigmaZetaPCA, mh4);
[Gvalue] = FMPM_Gval(NoSetup,Tempmat, InvSigmats, GG);

LowerBd=zeros(p,L0,2);
UpperBd=zeros(p,L0,2);
alphas=[0.05 0.01];

for pii=1:p;
Gkl=squeeze(Gvalue(pii,:,:));
Cklmat=zeros(1,GG);
for gii=1:GG
Cklmat(gii)=max(abs(squeeze(Gkl(:,gii))));
end
for ii=1:2
 Ckl=quantile(Cklmat,1-alphas(ii)); 
LowerBd(pii,:,ii)=squeeze(efitBetas(pii,:))-Ckl/sqrt(n)*ones(1,L0);
UpperBd(pii,:,ii)=squeeze(efitBetas(pii,:))+Ckl/sqrt(n)*ones(1,L0);
end
end
