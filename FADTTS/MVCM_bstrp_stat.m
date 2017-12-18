function [Gstat,Lstat] = MVCM_bstrp_stat(arclength,Xdesign,efitBetas,eSigEta,Cdesign,B0vector) 
% 
% MVCM_bstrp_stat is to implement Zhu's (2010) method of calculating the test statitics and p value using generated sample in MVCM
%
% Input:
%     arclength     - a L0 x 1 col vector of the arclength from one end to the other end
%     Xdesign       - a n x p normalized design matrix
%     efitBetas  - a p x L0 x m matrix of the estimated coefficients based on simulated data
%     Cdesign       - a r x mp matrix for characterizing the r linear constriants among mp parameters 
%     eSigEta    - a m x m x L0 x L0 matrix of the covariance of etas based on simulated data
%     B0vector      - a r x L0 vector for hypothesis testing
% Output:
%     Lstat      - a L0 x 1 columan vector of test statistics for each location based on simulated data
%     Gstat      - an over all test statitsics based on simulated data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%     [ResEta,efitEta,efitEta0,SigEta,h] = MVCM_sif(arclength,ResYdesign,m,kstr) 
% before you use MVCM_bstrp_stat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% April 4, 2010 @ AA
%     

L0=length(arclength);
n=size(Xdesign,1);
SimefitBetas=efitBetas;
SimeSigEta=eSigEta;
Omegax=pinv(Xdesign'*Xdesign);
SimLstat=zeros(L0,1);

for L0ii=1:L0
    TempBeta0=squeeze(SimefitBetas(:,L0ii,:));
    dd=Cdesign*TempBeta0(:)-B0vector(:,L0ii);
    SimLstat(L0ii)=dd'*pinv(Cdesign*(kron(squeeze(SimeSigEta(:,:,L0ii,L0ii)),Omegax))*Cdesign')*dd;
end
 
SimGstat=sum(SimLstat);
%SimGstat=(arclength(2:L0)-arclength(1:(L0-1)))'*(SimLstat(1:(L0-1))+SimLstat(2:L0))/2;
Gstat=SimGstat;
Lstat=SimLstat;
end