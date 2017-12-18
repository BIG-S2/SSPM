function [SimYdesign] = FMPM_grs(NoSetup, Indicator,efitBetas,efitEtas,ResEtas,Xdesign)
% 
% FMPM_grs is to generating samples
%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements
%     Indicator  - nx1 cell, where  Indicator{i}=k1:k2 is a r_ix1
%                  vector, contains the index for each subject i 
%     efitBetas  - a p x L0 matrix of estimated coefficients
%     efitEtas   - a sumri x L0  matrix of estimated etas
%     ResEtas    - a sumri x L0  matrix of epsilons
%     Xdesign    - a sumri x p matrix of covariates  
% Output:
%     SimYdesign      - a sumri x L0 matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude         
n=NoSetup(1);
L0=NoSetup(2);

taus=cell(n,1);
taus1=cell(n,1);

for i=1:n
    tau=normrnd(0,1);
    taus{i}=tau*ones(size(Indicator{i},1),L0); 
    tau1=normrnd(0,1,size(Indicator{i},1),L0);
    taus1{i}=tau1;
end

SimYdesign=Xdesign*efitBetas+cell2mat(taus).*efitEtas+cell2mat(taus1).*ResEtas;
  