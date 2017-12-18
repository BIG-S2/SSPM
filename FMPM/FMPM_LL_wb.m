function [efitBetas,efitBetas1,InvSigmats,efitYdesign, ResYdesign] = FMPM_LL_wb(NoSetup,arclength,Xdesign,Ydesign,mh,kstr)
% 
% FMPM_LL_wb is to implemen method of local polynomial kernel smoothing (order = 1) with preselected bandwidth
%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements (multivariate responses)
%     arclength    - L0x1 col vector of the arclength from one end to the other end
%     Xdesign    - a (\sum_i=1^n ri) x p normalized design matrix
%     Ydesign    - a ((\sum_i=1^n ri)xL0)  matrix of fiber bundle diffusion properties 
%     mh         - a 1 x 1 vector of optimal bandwidth
%     kstr       - Kernel function  
% Output:
%     efitBetas    - a p x L0  matrix of estimated coefficients
%     efitBetas1   - a 2*p x L0  matrix of estimated coefficients
%     InvSigmats   - a p*2 x p*2 x L0  matrix
%     efitYdesign  - a (\sum_i=1^n ri) x L0  matrix of estimated curves 
%     ResYdesign  - a (\sum_i=1^n ri) x L0  matrix of estimated curves 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright  by Y. Yuan 2014 at St. Jude  
    

if nargin<6
  kstr='exp(-.5*t.^2)';
end   

L0=NoSetup(2);
p=size(Xdesign,2);
sumri=NoSetup(4);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

InvSigmats=zeros(2*p,2*p,L0);
efitBetas=zeros(p,L0);
efitYdesign=zeros(sumri,L0);
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




