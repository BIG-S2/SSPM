function [SimYdesign] = MVCM_grs(efitBetas,efitEtas,ResEtas,Xdesign)
% 
% MVCM_grs is to implement Zhu's (2010) method of generating samples in MVCM
%
% Input:
%     efitBetas  - a p x L0 x m matrix of estimated coefficients
%     efitEtas     - a n x L0 x m matrix of estimated etas
%     ResEtas      - a n x L0 x m matrix of epsilons
%     Xdesign      - a n x p matrix of covariates  
% Output:
%     Ydesign      - a m x L_0 x n matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% April 3, 2010 @ AA
%            
[n L0 m]=size(ResEtas);
SimYdesign=zeros(n,L0,m);

Temptaus=normrnd(0,1,n,L0+1);
taus=Temptaus(:,1)*ones(1,L0);
taumat=Temptaus(:,2:end);
for mii=1:m
    SimYdesign(:,:,mii)=Xdesign*squeeze(efitBetas(:,:,mii))+taus.*squeeze(efitEtas(:,:,mii))+taumat.*squeeze(ResEtas(:,:,mii));
 end
   
end 
 
 

 


  