function [efitBetas,efitBetas1]= FMPM_LL_Ref_wb(NoSetup,Indicator, arclength,Xdesign,Ydesign, Zdesign,Sigma,SigmaEta, mh,kstr)
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
%                   properties
%     Zdesign    - (\sum_i=1^n ri) x p_z random normalized design matrix
%     Sigma      -  p_zxL0x pzxL0 of optimal covariance function 
%     SigmaEta   -  L0x1 covariance matrix for Eta
%     mh         -  bandwidth
%     kstr       -  Kernel function  
%%
% Output:
%     efitBetas     - p x L0  matrix of estimated coefficients
%     efitBetas1    - 2*p x L0  matrix of estimated coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude

if nargin<10, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
sumri=NoSetup(4);
pz=NoSetup(5);

vYdesign=zeros(sumri,L0);
tempv=zeros(p,p,L0,n);
     for nii=1:n
         ri=size(Indicator{nii},1);
         X1=Zdesign(Indicator{nii},:);
         X11=Xdesign(Indicator{nii},:);
         
         tempDD=zeros(ri^2,L0);
         
         for L0ii=1:L0
        Xs= X1*reshape(Sigma(:,L0ii),pz,pz)* X1'+SigmaEta(L0ii)*eye(ri);  

         tempDD(:,L0ii)=reshape(pinv(Xs), ri^2,1);

         vYdesign(Indicator{nii},L0ii)= pinv(Xs)*Ydesign(Indicator{nii},L0ii);  
     
         tempv(:,:,L0ii,nii)  =X11'*pinv(Xs)*X11;  
          end
          
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
