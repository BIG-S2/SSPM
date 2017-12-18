function [Lstat,Gstat] = FMPM_bstrp_stat(arclength,efitBetas,InvXCovX,Cdesign,B0vector) 
% 
% FMPM_bstrp_stat is to  calculate the test statitics and p value
%
% Input:
%     arclength     - a L0 x 1 col vector of the arclength from one end to 
%                     the other end
%     efitBetas     - a p x L0 matrix of the estimated coefficients based 
%                     on simulated data
%     Cdesign       - a r x p matrix for characterizing the r linear 
%                     constriants among p parameters 
%     InvXCovX      - p x p x L0 pinv(XCovX')
%     B0vector      - a r x L0 vector for hypothesis testing
% Output:
%     Lstat         - a L0 x 1 columan vector of test statistics for each 
%                     location based on simulated data
%     Gstat         - an over all test statitsics based on simulated data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude  


L0=length(arclength);
Lstat=zeros(L0,1);

for L0ii=1:L0
    TempBeta0=efitBetas(:,L0ii);
    dd=Cdesign*TempBeta0-B0vector(:,L0ii);
    Lstat(L0ii)=dd'*pinv(Cdesign*InvXCovX(:,:,L0ii)*Cdesign')*dd;  
end
 
Gstat=sum(Lstat);
end