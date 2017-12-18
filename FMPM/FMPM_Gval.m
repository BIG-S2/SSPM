function [Gvalue] = FMPM_Gval(NoSetup,Tempmat, InvSigmats, GG)
% 
% Input:
%   NoSetup    - col vector of [n, L0, p, sumri,p_z] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements (multivariate responses)
%                  p_z = number of covariates correspoinding to random
%                  effects
%  InvSigmats    2*px2*p*L0
%  Tempmat      - nx2pxL0
%     GG           - times of simulation 
% Output:
%     Gvalue       - a p x L0 x GG matrix of simulated matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Copyright  by Y. Yuan 2014 at St. Jude    


n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);

   
 Gvalue=zeros(p,L0,GG); 
for gii=1:GG
    taumat=randn(1,n); 
    for ii=1:L0
        KXRZmat=(taumat*Tempmat(:,:,ii))';
        eGval=squeeze(InvSigmats(:,:,ii))*KXRZmat;
        temp=sqrt(n)*eGval;
        Gvalue(:,ii,gii)=temp(1:p);
    end   
    
end

