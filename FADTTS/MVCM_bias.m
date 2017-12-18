function [ebiasBetas] = MVCM_bias(NoSetup,arclength,Xdesign,Ydesign,InvSigmats,mh,kstr) 
% 
% MVCM_ht_stat is to implement Zhu's (2010) method of calculating the estimate bias in MVCM
%
% Input:
%     NoSetup      - col vector of [n, L0, p, m] 
%                    n = sample size 
%                    L0 = number of location of fiber tract 
%                    p = number of covariates 
%                    m = number of features 
%     arclength     - a L0 x 1 col vector of the arclength from one end to the other end
%     Xdesign       - a n x p normalized design matrix
%     Ydesign       - a m x L0 x n matrix of fiber bundle diffusion properties 
%     InvSigmats    - a 2p x 2p x L0 x m matrix
%     mh            - a 1 x m column vector of pilot optimal bandwidths
% Output:
%     ebiasBetas    - a p x L0 x m matrix of the bias of the estimated betas
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%     [ResEta,efitEta,efitEta0,SigEta,h] = MVCM_sif(arclength,ResYdesign,m,kstr) 
% before you use MVCM_bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 29, 2010 @ AA
%     

if nargin<7, 
  kstr='exp(-.5*t.^2)';
end 

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
m=NoSetup(4);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

[efitBetas3] = MVCM_lpks_wb2(NoSetup,arclength,Xdesign,Ydesign,mh);

ebiasBetas=zeros(p,L0,m);

for mii=1:m
    h=mh(mii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    Tempmat2=zeros(L0,2*p,L0);
    Tempmat3=zeros(L0,2*p,L0);
    efitBeta2=kron(eye(p),[0 0 1 0])*efitBetas3(:,:,mii);
    efitBeta3=kron(eye(p),[0 0 0 1])*efitBetas3(:,:,mii);
    
    TmatA=repmat(Tmat,[1,1,p]);
    for nii=1:n
        for L0ii=1:L0
            TempA=repmat(Xdesign(nii,:),[1,L0]);
            TempB=kron(Tmat(:,L0ii)',ones(1,p));
            TempC=[TempA;TempA.*TempB];
            TempD=reshape(TempC(:),[2*p,L0]);
            TempDD=permute(TempD,[2 1]);
            AA=Kmat(:,L0ii)'.*(Xdesign(nii,:)*(efitBeta2.*(squeeze(TmatA(:,L0ii,:))').^2+efitBeta3.*(squeeze(TmatA(:,L0ii,:))').^3));
            Tempmat2(:,:,L0ii)=TempDD.*repmat(AA',[1,2*p]);
           
        end
        Tempmat3=Tempmat3+Tempmat2;
    end

    KXBZmat=squeeze(sum(Tempmat3));

    InvSigmat=squeeze(InvSigmats(:,:,:,mii));
    ebiasBeta0=zeros(2*p,L0);
    for L0ii=1:L0
        ebiasBeta0(:,L0ii)=squeeze(InvSigmat(:,:,L0ii))*KXBZmat(:,L0ii);
    end
    
    ebiasBetas(:,:,mii)=kron(eye(p),[1 0])*ebiasBeta0;
end 

end