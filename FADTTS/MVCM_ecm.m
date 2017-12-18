function [eSigE] = MVCM_ecm(arclength,ResEtas,kstr) 
% 
% MVCM_ecm is to implement Zhu's (2010) method of smoothing individual function without preselected bandwidth in MVCM
%
% Input:
%     arclength  - col vector of the arclength from one end to the other end
%     ResEtas    - a n x L0 x m matrix of residuals
%     kstr       - Kernel function  
% Output:
%     eSigE      - a L0 x m x m matrix of the covariance of residuals 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%    function [ResEta,efitEta,efitEta0,SigEta,h] = MVCM_sif(arclength,ResYdesign,kstr) 
% before you use MVCM_ecm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 28, 2010 @ AA
%     

if nargin<3, 
  kstr='exp(-.5*t.^2)';
end 

[n L0 m]=size(ResEtas);
xrange=max(arclength)-min(arclength);
nh=floor(max(30, L0/2));
hmin=1.0*xrange/L0; 
hmax=xrange/8 ; 
vh=logspace(log10(hmin),log10(hmax),nh);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';
GCVs=zeros(nh,1);

Emat=zeros(n,L0,m,m);
for nii=1:n
    for L0ii=1:L0
        Emat(nii,L0ii,:,:)=squeeze(ResEtas(nii,L0ii,:))*(squeeze(ResEtas(nii,L0ii,:)))';
    end
end

tSigE=reshape(sum(Emat(:,:)),[L0,m,m]);
InvtSigE=tSigE*0.0;
for L0ii=1:L0
    InvtSigE(L0ii,:,:)=pinv(squeeze(tSigE(L0ii,:,:)))/(n-m);
end

for nhii=1:nh
    h=vh(nhii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    sumKmat=sum(Kmat);
    meanKmat=Kmat./(sumKmat'*ones(1,L0));
    
    CC=zeros(n,m,m,L0);
    for nii=1:n
        for L0ii=1:L0
            AA=squeeze(Emat(nii,:,:,:));
            BB=repmat(meanKmat(:,L0ii),[1 m m]);
            TempCC=sum(AA(:,:).*BB(:,:));
            CC(nii,:,:,L0ii)=reshape(TempCC,[m,m]);
        end
    end
 
    BB=squeeze(sum(CC,1));
    GCVmat=zeros(n,L0);
    for nii=1:n
        jSigE=BB-squeeze(CC(nii,:,:,:));
        for L0ii=1:L0
            GCVmat(nii,L0ii)=(trace(squeeze(Emat(nii,L0ii,:,:))-squeeze(jSigE(:,:,L0ii))*squeeze(InvtSigE(L0ii,:,:))))^2;
        end
    end
    GCVs(nhii)=sum(sum(GCVmat))/n/L0;
end

[MinGCV,flag]=min(GCVs);
h=vh(flag);

Tmat=Tmat0/h;
t=Tmat;
Kmat=eval(kstr)/h;
sumKmat=sum(Kmat);
meanKmat=Kmat./(sumKmat'*ones(1,L0));

eSigE=zeros(L0,m,m);
for L0ii=1:L0
    TempBB=zeros(m,m);
    for nii=1:n
        TempCC=zeros(m,m);
        for L0jj=1:L0
            TempCC=TempCC+squeeze(Emat(nii,L0jj,:,:))*meanKmat(L0jj,L0ii);
        end
        TempBB=TempBB+TempCC;
    end
    eSigE(L0ii,:,:)=TempBB/(n-m);
end
   
end