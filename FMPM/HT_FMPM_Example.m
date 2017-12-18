clear all;
close all;
fclose all;
addpath('/');
workDir='/RealData_NDARTBSS';
PREfile='RealData_NDAR';
APENfile='mat';
DataFile=sprintf('%s/%s.%s', workDir,PREfile,APENfile);
load(DataFile); 
Xdesign=[Xdesign(:,1:2), log(Xdesign(:,3)), log(Xdesign(:,3)).^2];
Zdesign= Xdesign(:,[1,3]);
[NoSetup, Xdesign, Zdesign]=FMPM_ReadNormalize(Xdesign,Zdesign,Indicator, arclength);
n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
sumri=NoSetup(4);
pz=NoSetup(5);
Cdesign=[0,0,1,0;0 0 0 1];
B0vector=zeros(2,L0);
ExpVar=0.99;
rng('default');
[mh1,~,~] = FMPM_LL_wob(NoSetup,Indicator,arclength,Xdesign,Ydesign);
[~,~,InvSigmats,~,~,ResYdesign,H0ResYdesign] = FMPM_LL_wb_H0(NoSetup,arclength,Xdesign,Ydesign,mh1,Cdesign,B0vector);
[SigmaZeta,SigmaZetad,SigmaXi,SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,ResYdesign);
[H0SigmaZeta,H0SigmaZetad,H0SigmaXi,H0SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,H0ResYdesign);
[mh0] = BiSmooth_LL_bs(NoSetup,arclength,SigmaXi,pz,0);
[SigmaXiP,SigmaXiE] = BiSmooth_LL(NoSetup,arclength,SigmaXi,SigmaXid,pz,mh0,0);
[mh00] = BiSmooth_LL_bs(NoSetup,arclength,SigmaZeta,1,0);
[SigmaZetaP,SigmaZetaE] = BiSmooth_LL(NoSetup,arclength,SigmaZeta, SigmaZetad,1,mh00,0);
[mh0] = BiSmooth_LL_bs(NoSetup,arclength,H0SigmaXi,pz,0);
[H0SigmaXiP,H0SigmaXiE] = BiSmooth_LL(NoSetup,arclength,H0SigmaXi,H0SigmaXid,pz,mh0,0);
[mh00] = BiSmooth_LL_bs(NoSetup,arclength,H0SigmaZeta,1,0);
[H0SigmaZetaP,H0SigmaZetaE] = BiSmooth_LL(NoSetup,arclength,H0SigmaZeta, H0SigmaZetad,1,mh00,0);
SigmaXiE=kron(eye(L0),SigmaXiE);
SigmaZetaE=SigmaZetaE*eye(L0);
H0SigmaXiE=kron(eye(L0),H0SigmaXiE);
H0SigmaZetaE=H0SigmaZetaE*eye(L0);

[~,~,~,~,~,~,SigmaXiPCA, SigmaZetaPCA] = FMPM_FPCA_Cov(NoSetup,SigmaXiP,SigmaZetaP,SigmaXiE,SigmaZetaE,ExpVar);
[~,~,~,~,~,~,H0SigmaXiPCA, H0SigmaZetaPCA] = FMPM_FPCA_Cov(NoSetup,H0SigmaXiP,H0SigmaZetaP,H0SigmaXiE,H0SigmaZetaE,ExpVar);
[mh4,GCVs,vh,vYdesign,vYdesignH0, tempv, tempvH0] = FMPM_LL_Ref_wob_H0(NoSetup,Indicator,arclength,Xdesign,Zdesign,Ydesign, SigmaXiPCA,SigmaZetaPCA, H0SigmaXiPCA,H0SigmaZetaPCA);
mh4=mh4/6;
[efitBetas,efitBetas1,efitBetasH0,efitBetas1H0,efitYdesign,efitYdesignH0, ResYdesign, H0ResYdesign, InvXCovX]= FMPM_LL_Ref_wb_H0(NoSetup,Indicator, arclength,Xdesign,Ydesign,vYdesign,tempv,mh4,Cdesign,B0vector);
efitBetas_2=efitBetas;

[Lstat,Gstat] = FMPM_bstrp_stat(arclength,efitBetas,InvXCovX,Cdesign,B0vector);
[ResEtasH0,efitEtasH0] = FMPM_LL_If_wob(NoSetup,arclength,H0ResYdesign);
GG=1000;
Smat=zeros(L0,GG);
SAmat=zeros(GG,1);
Pvalue=zeros(L0,1);
PvalueG=0;
for gii=1:GG
[gYdesign] = FMPM_grs(NoSetup,Indicator,efitBetasH0,efitEtasH0,ResEtasH0,Xdesign);
[gefitBetas,gefitBetas1]= FMPM_LL_Ref_wb(NoSetup,Indicator,arclength,Xdesign,gYdesign,Zdesign,SigmaXiPCA,SigmaZetaPCA,mh4);
[gLstat,gGstat] = FMPM_bstrp_stat(arclength,gefitBetas,InvXCovX,Cdesign,B0vector);
Smat(:,gii)=gLstat;
SAmat(gii)=gGstat;
MaxSn=max(gLstat);
Pvalue=Pvalue+( (MaxSn*ones(L0,1))>Lstat);
PvalueG=PvalueG+(gGstat>Gstat);
end
PvalL=Pvalue/GG;
PvalG=PvalueG/GG;


