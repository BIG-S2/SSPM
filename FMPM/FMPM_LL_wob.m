function [mh,GCVs,vh] = FMPM_LL_wob(NoSetup,Indicator, arclength,Xdesign,Ydesign,kstr)
% 
% FMPM_LL_wob is to implement local polynomial kernel smoothing of varying 
% coefficient function  without preselected bandwidth
%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements=\sum_i=1^n r_i
%     Indicator  - nx1 cell, where  Indicator{i}=k1:k2 is a r_ix1
%                  vector, contains the index for each subject i in array
%     arclength  - L0x1 col vector of the arclength from one end to the 
%                  other end
%     Xdesign    - (\sum_i=1^n r_i) x p normalized design matrix
%     Ydesign    - (\sum_i=1^n r_i)xL0  matrix of fiber bundle diffusion 
%                  properties 
%     kstr       - Kernel function  
% Output:
%     mh         - optimal bandwidth
%     GCVs       - nh x1 matrix of GCVs
%     vh         - 1 x nh vector of bandwidth
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    Copyright  by Y. Yuan 2014 at St. Jude
    
if nargin<6, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=size(Xdesign,2);
sumri=NoSetup(4);

xrange=max(arclength)-min(arclength);
nh=floor(L0/2);
hmin=1.0*xrange/L0; 
hmax=xrange/8; 
efitYdesign=Ydesign*0;
GCVs=zeros(nh,1);
vh=logspace(log10(hmin),log10(hmax),nh);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

for nhii=1:nh
    h=vh(nhii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
  
Tempmat0=zeros(2*p,2*p,L0,sumri);
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
            Tempmat0(:,:,L0ii,nii)=TempEE'*TempDD;
        end
    end   
    
        AA=repmat(Ydesign,[1 1 2*p L0]);
        Tempmat10=Tempmat(:,:).*AA(:,:); 
        Tempmat11=reshape(Tempmat10,[sumri,L0,2*p,L0]);        
        Tempmat1=squeeze(sum(Tempmat11,2));        
        
        for L0ii=1:L0
            Tempmat2=sum(squeeze(Tempmat0(:,:,L0ii,:)),3);
            Tempmat3=sum(squeeze(Tempmat1(:,:,L0ii)),1);
      
            for nii=1:n 
                Tempmat22=Tempmat2-sum(squeeze(Tempmat0(:,:,L0ii,Indicator{nii})),3);
                Tempmat33=Tempmat3-sum(squeeze(Tempmat1(Indicator{nii},:,L0ii)),1);
                AA1=pinv(Tempmat22)*Tempmat33';
                efitYdesign(Indicator{nii},L0ii)= Xdesign(Indicator{nii},:)*kron(eye(p),[1 0])*AA1;
            end
        end
        GCVs(nhii)=sum(sum((Ydesign-efitYdesign).^2))/(sumri*L0); 
end

[~,flag]=min(GCVs);
mh=vh(flag);