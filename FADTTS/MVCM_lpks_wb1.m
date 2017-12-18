function [efitBetas,efitBetas1,InvSigmats,efitYdesign] = MVCM_lpks_wb1(NoSetup,arclength,Xdesign,Ydesign,mh,kstr)
% 
% MVCM_lpks_wb1 is to implement Zhu's (2010) method of local polynomial kernel smoothing (order = 1) with preselected bandwidth in MVCM
%
% Input:
%     NoSetup      - col vector of [n, L0, p, m] 
%                    n = sample size 
%                    L0 = number of location of fiber tract 
%                    p = number of covariates 
%                    m = number of features 
%     arclength    - col vector of the arclength from one end to the other end
%     Xdesign      - a n x p normalized design matrix. 
%     Ydesign      - a n x L0 x m matrix of fiber bundle diffusion properties 
%     mh           - a 1 x m vector of optimal bandwidth
%     kstr         - Kernel function  
% Output:
%     efitBetas    - a p x L0 x m matrix of estimated coefficients
%     efitBetas1   - a 2*p x L0 x m matrix of estimated coefficients
%     InvSigmats   - a p*2 x p*2 x L0 x m matrix
%     efitYdesign  - a n x L0 x m matrix of estimated curves 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%    function [NoSetup, arclength, Xdesign,Ydesign]=MVCM_read(tractdata, designdata,  diffusionFiles, nofeatures, featuresname)  
% before you use MVCM_lpks_wb1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 27, 2010 @ AA
%     

if nargin<6, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
m=NoSetup(4);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

InvSigmats=zeros(2*p,2*p,L0,m);
efitBetas=zeros(p,L0,m);
efitYdesign=zeros(n,L0,m);
efitBetas1=zeros(p*2,L0,m);

for mii=1:m
    
    h=mh(mii); 
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;

    Sigmat=zeros(2*p,2*p,L0);
    Tempmat0=zeros(2*p,2*p,L0);
    Tempmat=zeros(n,L0,2*p,L0);
    for nii=1:n
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
        InvSigmats(:,:,L0ii,mii)=pinv(Sigmat(:,:,L0ii));
    end
      
    AA=repmat(squeeze(Ydesign(:,:,mii)),[1 1 2*p]);
    for L0ii=1:L0
        BB=squeeze(Tempmat(:,:,:,L0ii));
        CC=AA(:,:).*BB(:,:);
        CCC=sum(CC);
        KXYZmat=sum(reshape(CCC,[L0,2*p]));
        efitBetas1(:,L0ii,mii)=squeeze(InvSigmats(:,:,L0ii,mii))*KXYZmat';
    end

    efitBeta=kron(eye(p),[1 0])*efitBetas1(:,:,mii);
    efitBetas(:,:,mii)=efitBeta;
    efitYdesign(:,:,mii)=Xdesign*efitBeta;
end

end

