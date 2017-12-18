function [efitBetas,efitBetas1,efitBetasH0,efitBetas1H0,efitYdesign,efitYdesignH0, ResYdesign, H0ResYdesign,InvXCovX]= FMPM_LL_Ref_wb_H0(NoSetup,Indicator, arclength,Xdesign,Ydesign,vYdesign,tempv, mh,Cdesign,B0vector,kstr)
%%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri,p_z] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements 
%                  p_z = number of covariates correspoinding to random
%                  effects
%     Indicator  - nx1 cell, where  Indicator{i}=k1:k2 is a r_ix1
%                  vector, contains the index for each subject i
%     arclength  - L0x1 col vector of the arclength from one end to the 
%                  other end
%     Xdesign    - (\sum_i=1^n ri) x p normalized design matrix
%     Ydesign    - ((\sum_i=1^n ri)xL0)  matrix of fiber bundle diffusion
%                  properties 
%     vYdesign   - sumri x L0 matrix
%     tempv      - (p,p,L0,n) arrary
%     Cdesign    - r x p matrix for characterizing the r linear constriants
%                  among p parameters 
%     B0vector   - r x L0 vector
%     mh         - bandwidth
%     kstr       - Kernel function  
%%
% Output:
%     efitBetas     - p x L0  matrix of estimated coefficients
%     efitBetas1    - 2*p x L0  matrix of estimated coefficients
%     efitBetasH0   - p x L0  matrix of estimated coefficients
%     efitBetas1H0  - 2*p x L0  matrix of estimated coefficients
%     efitYdesign   - (\sum_i=1^n ri) x L0  matrix of estimated curves 
%     efitYdesignH0 - (\sum_i=1^n ri) x L0  matrix of estimated curves 
%     ResYdesign    - (\sum_i=1^n ri) x L0  matrix of estimated curves
%     H0ResYdesign  - (\sum_i=1^n ri) x L0  matrix of estimated curves under H0 
%     InvXCovX      - p x p x L0 pinv(XCovX')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude  

if nargin<11, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);

 InvXCovX=zeros(p,p,L0);
 for L0ii=1:L0
     InvXCovX(:,:,L0ii)=pinv(sum(squeeze(tempv(:,:,L0ii,:)), 3));
 end
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';



    h=mh;
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
  
Tempmat0=zeros(2*p*2*p,L0,L0,n);
    Tempmat=zeros(n,2*p,L0);
  
    for nii=1:n
       for L0ii=1:L0
       Temp3= Xdesign(Indicator{nii},:)';
       Temp4=[Kmat(:,L0ii)';Kmat(:,L0ii)'.*Tmat(:,L0ii)'];
       Temp5=kron(Temp4, Temp3);   
       
       Temp51=vYdesign(Indicator{nii},:);  
       Temp52=Temp51(:); 
       Tempmat(nii,:,L0ii)=Temp5*Temp52;

       Temp1=repmat(Kmat(:,L0ii)',4,1).*[ones(1,L0);Tmat(:,L0ii)';Tmat(:,L0ii)';Tmat(:,L0ii)'.*Tmat(:,L0ii)']; 
       Temp2=repmat(Temp1(:), 1,p^2)'; 
       
       Temp6=repmat(reshape(squeeze(tempv(:,:,:,nii)), p^2,L0)', 1,4)';  
       Temp7=reshape(Temp6, p^2, 4*L0);
       Tempmat0(:,:,L0ii,nii)=reshape(Temp2.*Temp7,4*p^2,L0);   
       end
    end

efitBetas1=zeros(2*p,L0);
InvSigmats=zeros(2*p,2*p,L0);
        for L0ii=1:L0
            Tempmat2=sum(sum(squeeze(Tempmat0(:,:,L0ii,:)),3),2);
            Tempmat3=sum(squeeze(Tempmat(:,:,L0ii)),1);
          
             tttt=reshape(Tempmat2,p,p,4);
             tttt1=[tttt(:,:,1), tttt(:,:,2);  tttt(:,:,3), tttt(:,:,4)]; 
             efitBetas1(:,L0ii)=pinv(tttt1)*Tempmat3';
             InvSigmats(:,:,L0ii)=pinv(tttt1); 
        end
        
efitBetas=efitBetas1(1:p,:);
efitYdesign=Xdesign*efitBetas;
ResYdesign=Ydesign-efitYdesign;  
efitBetas1H0=efitBetas1*0.0;

AA=[eye(p),zeros(p,p)];
CdesignNew=Cdesign*AA;

for L0ii=1:L0
    VecBeta10=efitBetas1(:,L0ii);
    InvSigs=InvSigmats(:,:,L0ii);
    efitBetas1H0(:,L0ii)=VecBeta10-InvSigs*CdesignNew'*pinv(CdesignNew*InvSigs*CdesignNew')*(CdesignNew*VecBeta10-B0vector(:,L0ii));
end

efitBetasH0=AA*efitBetas1H0;

efitYdesignH0=Xdesign*efitBetasH0;

H0ResYdesign=Ydesign-efitYdesignH0;

