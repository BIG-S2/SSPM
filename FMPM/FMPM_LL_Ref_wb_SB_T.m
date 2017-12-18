function [efitBetas,efitBetas1,efitYdesign, ResYdesign,InvSigmats, Tempmat]= FMPM_LL_Ref_wb_SB_T(NoSetup,Indicator, arclength,Xdesign,Zdesign,Ydesign,vYdesign,tempv,Sigma,SigmaEta, mh,kstr)
%%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri,p_z] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements (multivariate responses)
%                  p_z = number of covariates correspoinding to random
%                  effects
%     Indicator  - nx1 cell, where  Indicator{i}=k1:k2 is a r_ix1
%                  vector, contains the index for each subject i in array
%                  Ydesign.r_i = number of time points for subject i.
%                  -----very important
%     arclength    - L0x1 col vector of the arclength from one end to the other end
%     Xdesign    - (\sum_i=1^n ri) x p normalized design matrix
%     Zdesign    - (\sum_i=1^n ri) x p_z random normalized design matrix
%     Ydesign    - ((\sum_i=1^n ri)xL0)  matrix of fiber bundle diffusion properties 
%     vYdesign   - a sumri x L0 matrix
%     tempv      - a (p,p,L0,n) arrary
%     Sigma     - a p_z^2xL0 of optimal covariance function 
%     SigmaEta  - a L0x1 covariance matrix for Eta
%     mh         - bandwidth
%     kstr       - Kernel function  
%%
% Output:
%     efitBetas    - a p x L0  matrix of estimated coefficients
%     efitBetas1   - a 2*p x L0  matrix of estimated coefficients
%     efitYdesign  - a (\sum_i=1^n ri) x L0  matrix of estimated curves 
%     ResYdesign   - a (\sum_i=1^n ri) x L0  matrix of estimated curves
%     InvSigmats   - 2*px2*p*L0
%     Tempmat      - nx2pxL0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude   
   

if nargin<12, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
sumri=NoSetup(4);
pz=NoSetup(5);

Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
            Tempmat2=sum(sum(squeeze(Tempmat0(:,:,L0ii,:)),3),2);%
            Tempmat3=sum(squeeze(Tempmat(:,:,L0ii)),1);
             tttt=reshape(Tempmat2,p,p,4);
             tttt1=[tttt(:,:,1), tttt(:,:,2);  tttt(:,:,3), tttt(:,:,4)]; 
             efitBetas1(:,L0ii)=pinv(tttt1)*Tempmat3';
             InvSigmats(:,:,L0ii)=pinv(tttt1); 
        end
        
efitBetas=efitBetas1(1:p,:);
efitYdesign=Xdesign*efitBetas;
ResYdesign=Ydesign-efitYdesign;













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vResYdesign=zeros(sumri,L0);
     for nii=1:n
         ri=size(Indicator{nii},1);
         X1=Zdesign(Indicator{nii},:);
         for L0ii=1:L0  
         Xs= X1*reshape(Sigma(:,L0ii),pz,pz)* X1'+SigmaEta(L0ii)*eye(ri); 
         vResYdesign(Indicator{nii},L0ii)= pinv(Xs)*ResYdesign(Indicator{nii},L0ii);  % sumrixL0
         end
     end

     
  for nii=1:n   
       for L0ii=1:L0
       Temp3= Xdesign(Indicator{nii},:)';
       Temp4=[Kmat(:,L0ii)';Kmat(:,L0ii)'.*Tmat(:,L0ii)']; %2*L0
       Temp5=kron(Temp4, Temp3);   %2px(ri*L0);
       Temp51=vResYdesign(Indicator{nii},:);  %ri*L0
       Temp52=Temp51(:); %ri*L0x1
       Tempmat(nii,:,L0ii)=Temp5*Temp52; % (n,2p,L0)
       end   
  end
    








