function [Gstat,Lstat] = MVCM_ht_stat(NoSetup,arclength,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas) 
% 
% MVCM_ht_stat is to implement Zhu's (2010) method of calculating the test statitics and p value in MVCM
%
% Input:
%     NoSetup      - col vector of [n, L0, p, m] 
%                    n = sample size 
%                    L0 = number of location of fiber tract 
%                    p = number of covariates 
%                    m = number of features 
%     arclength     - a L0 x 1 col vector of the arclength from one end to the other end
%     Xdesign       - a n x p normalized design matrix
%     efitBetas     - a p x L0 x m matrix of the estimated coefficients
%     Cdesign       - a r x mp matrix for characterizing the r linear constriants among mp parameters 
%     eSigEta       - a m x m x L0 x L0 matrix of the covariance of etas
%     B0vector      - a r x L0 vector for hypothesis testing
%     ebiasBetas    - a p x L0 x mmatrix of the bias of the estimated betas
% Output:
%     Lstat         - a L0 x 1 columan vector of test statistics for each location 
%     Gstat         - an over all test statitsics
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%     [ResEta,efitEta,efitEta0,SigEta,h] = MVCM_sif(arclength,ResYdesign,m,kstr) 
% before you use MVCM_ht_stat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 29, 2010 @ AA
%     

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
m=NoSetup(4);

deltaBeta=zeros(m*p,L0);

for mii=1:m
    deltaBeta((mii-1)*p+1:mii*p,:)=efitBetas(:,:,mii)-ebiasBetas(:,:,mii);        
end 

dd=Cdesign*deltaBeta-B0vector;
Omegax=pinv(Xdesign'*Xdesign);
Lstat=zeros(L0,1);
for L0ii=1:L0
    Lstat(L0ii)=dd(:,L0ii)'*pinv(Cdesign*(kron(squeeze(eSigEta(:,:,L0ii,L0ii)),Omegax))*Cdesign')*dd(:,L0ii);
end

Gstat=sum(Lstat);
%Gstat=(arclength(2:L0)-arclength(1:(L0-1)))'*(Lstat(1:(L0-1))+Lstat(2:L0))/2;

end