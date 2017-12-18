function [efitBetas,efitBetas1,InvSigmats,efitYdesign,H0efitBetas, ResYdesign, H0ResYdesign] = FMPM_LL_wb_H0(NoSetup,arclength,Xdesign,Ydesign,mh,Cdesign,B0vector,kstr)
% 
% FMPM_LL_wb_H0 is to implement method of local polynomial kernel smoothing 
% with selected bandwidth
%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements
%     arclength  - L0x1 col vector of the arclength from one end to the 
%                  other end
%     Xdesign    - (\sum_i=1^n ri) x p normalized design matrix
%     Ydesign    - ((\sum_i=1^n ri)xL0)  matrix of fiber bundle diffusion
%                  properties 
%     mh         - optimal bandwidth 
%     Cdesign    - r x p matrix for characterizing the r linear constriants
%                  among p parameters 
%     B0vector   - r x L0 vector
%     kstr       - Kernel function 
% Output:
%     efitBetas    - p x L0  matrix of estimated coefficients
%     H0efitBetas  - p x L0  matrix of estimated coefficients under H0
%     efitBetas1   - 2*p x L0  matrix of estimated coefficients
%     InvSigmats   - p*2 x p*2 x L0  matrix
%     efitYdesign  - (\sum_i=1^n ri) x L0  matrix of estimated curves 
%     ResYdesign   - (\sum_i=1^n ri) x L0  matrix of estimated curves 
%     H0ResYdesign - (\sum_i=1^n ri) x L0  matrix of estimated curves under
%                    H0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude    

if nargin<8
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
sumri=NoSetup(4);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

InvSigmats=zeros(2*p,2*p,L0);
efitBetas1=zeros(p*2,L0);
    
    Tmat=Tmat0/mh;
    t=Tmat;
    Kmat=eval(kstr)/mh;
    
  Sigmat=zeros(2*p,2*p,L0);
   Tempmat0=zeros(2*p,2*p,L0);
    Tempmat=zeros(sumri,L0,2*p,L0);
    for nii=1:sumri
        for L0ii=1:L0
            TempA=repmat(Xdesign(nii,:),[1,L0]);
            TempB=kron(Tmat(:,L0ii)',ones(1,p));
            TempC=[TempA;TempA.*TempB];
            TempD=reshape(TempC(:),[2*p,L0]);
            TempDD=permute(TempD,[2 1]);   
            TempEE=TempDD.*repmat(Kmat(:,L0ii),[1,2*p]);
            Tempmat(nii,:,:,L0ii)=TempEE;
            Tempmat0(:,:,L0ii)=TempEE'*TempDD;
        end
         Sigmat=Sigmat+Tempmat0;
    end
    

    for L0ii=1:L0
        InvSigmats(:,:,L0ii)=pinv(Sigmat(:,:,L0ii));
    end
      
    AA=repmat(Ydesign,[1 1 2*p]);
    for L0ii=1:L0
        BB=squeeze(Tempmat(:,:,:,L0ii));
        CC=AA(:,:).*BB(:,:);
        CCC=sum(CC);
        KXYZmat=sum(reshape(CCC,[L0,2*p]));
        efitBetas1(:,L0ii)=squeeze(InvSigmats(:,:,L0ii))*KXYZmat';
    end

    efitBeta=kron(eye(p),[1 0])*efitBetas1;
    efitBetas=efitBeta;
    efitYdesign=Xdesign*efitBeta;
    ResYdesign=Ydesign-efitYdesign;
    
H0efitBetas1=efitBetas1*0.0;

AA1=kron(eye(p),[1 0]);
CdesignNew=Cdesign*AA1;

for L0ii=1:L0
    VecBeta10=efitBetas1(:,L0ii);
    InvSigs=InvSigmats(:,:,L0ii);
    H0efitBetas1(:,L0ii)=VecBeta10-InvSigs*CdesignNew'*pinv(CdesignNew*InvSigs*CdesignNew')*(CdesignNew*VecBeta10-B0vector(:,L0ii));
end

H0efitBetas=AA1*H0efitBetas1;
H0efitYdesign=Xdesign*H0efitBetas; 
H0ResYdesign=Ydesign-H0efitYdesign;