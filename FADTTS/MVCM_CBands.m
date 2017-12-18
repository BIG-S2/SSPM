function [CBands] = MVCM_CBands(n,alpha,Gvalue,efitBetas,ebiasBetas) 
% 
% MVCM_cb_Gval is to implement Zhu's (2010) method of calculating confidence bands in MVCM
%
% Input:
%     Gvalue       - a m x p x L0 x GG matrix of simulated matrix
%     efitBetas    - a p x L0 x m matrix of estimated coefficients
%     ebiasBetas   - a p x L0 x m matrix of the bias of the estimated betas
% Output:
%     CBands       - a 2p x L0 x m matrix of estimated confidence bands
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%           [Gvalue] = MVCM_cb_Gval(arclength,Xdesign,ResYdesign,InvSigmats,mh,GG,kstr) 
% before you use MVCM_CBands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% September 3, 2010 @ AA
%     
        
        p=size(efitBetas,1);
        L0=size(efitBetas,2);
        m=size(efitBetas,3);
        GG=size(Gvalue,4);
        CBands=zeros(2*p,L0,m);
        for mii=1:m
            for pii=1:p
                Gkl=squeeze(Gvalue(mii,pii,:,:));
                Cklmat=zeros(1,GG);
                for gii=1:GG
                    Cklmat(gii)=max(abs(squeeze(Gkl(:,gii))));
                end
                Ckl=quantile(Cklmat,1-alpha);
                CBands(2*pii-1,:,mii)=squeeze(efitBetas(pii,:,mii)-ebiasBetas(pii,:,mii))-Ckl/sqrt(n)*ones(1,L0);
                CBands(2*pii,:,mii)=squeeze(efitBetas(pii,:,mii)-ebiasBetas(pii,:,mii))+Ckl/sqrt(n)*ones(1,L0);
            end
        end
end