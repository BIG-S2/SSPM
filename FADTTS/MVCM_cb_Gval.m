function [Gvalue] = MVCM_cb_Gval(arclength,Xdesign,ResYdesign,InvSigmats,mh,GG,kstr) 
% 
% MVCM_cb_Gval is to implement Zhu's (2010) method of approximating G value in MVCM
%
% Input:
%     arclength    - col vector of the arclength from one end to the other end
%     Xdesign      - a n x p matrix of covariates
%     ResYdesign   - a n x L0 x m matrix of difference between fiber bundle diffusion properties and fitted fiber bundle diffusion properties
%     InvSigmats   - a 2p x 2p x L0 x m matrix
%     mh           - a 1 x m column vector of pilot optimal bandwidths
%     kstr         - Kernel function 
%     GG           - times of simulation 
% Output:
%     Gvalue       - a m x p x L0 x GG matrix of simulated matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%     function [efitYdata,efitBeta,h] = MVCM_lpks_wob(NoSetup,arclength,Xdesign,Ydata,korder,kstr)
% before you use MVCM_sif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% April 5, 2010 @ AA
%     

if nargin<7, 
  kstr='exp(-.5*t.^2)';
end 

p=size(Xdesign,2);
[n L0 m]=size(ResYdesign);
Gvalue=zeros(m,p,L0,GG);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';
Tempmat=zeros(n,2*p,L0,m);

for mii=1:m
    h=mh(mii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
  
    for nii=1:n
        for L0ii=1:L0
            AA=zeros(2*p,1);
            for L0jj=1:L0
                Tempt=[Xdesign(nii,:);Xdesign(nii,:)*Tmat(L0jj,L0ii)];
                AA=AA+Tempt(:)*Kmat(L0jj,L0ii)*ResYdesign(nii,L0jj,mii);
            end
            Tempmat(nii,:,L0ii,mii)=AA;
        end
    end
end

for gii=1:GG    
    taumat=randn(1,n);
    for mii=1:m
        eGval=zeros(p*2,L0);
        for L0ii=1:L0
            KXRZmat=(taumat*Tempmat(:,:,L0ii,mii))';
            eGval(:,L0ii)=squeeze(InvSigmats(:,:,L0ii,mii))*KXRZmat;
        end
        Gvalue(mii,:,:,gii)=sqrt(n)*kron(eye(p),[1 0])*eGval;
    end    
end    


end