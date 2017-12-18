function [NXi,NZeta,eveXi,eveZeta,evaXi,evaZeta,sigmaXiPCA, SigmaZetaPCA]=FMPM_FPCA_Cov(NoSetup,SigmaXiP,SigmaZetaP,SigmaXiE,SigmaZetaE,ExpVar,NXi, NZeta)
% FMPM_FPCA_Cov functional principal component analysis in FMPM
%
% Input:
%     NoSetup      - col vector of [n, L0, p,sumri] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements 
%                  p_z = number of covariates correspoinding to random
%                  effects
%     Indicator  - nx1 cell, where  Indicator{i}=k1:k2 is a r_ix1
%                  vector, contains the index for each subject i 
%     SigmaXiP     - L0*p_z x L0*p_z matrix 
%     SigmaZetaP   - L0 x L0 matrix 
%     Sigma2       - L0 x1 vector
%     ExpVar       - Specified expained variance
%     Zdesign      - a (\sum_i=1^n ri) x p_z normalized design matrix
%     ResYdesign   - a ((\sum_i=1^n ri)xL0)   matrix 
%     NXi          - number of eigenfunctions with respect to SigmaXiP 
%                    retained for the analysis
%     NZeta        - number of eigenfunctions with respect to SigmaZetaP 
%                    retained for the analysis
% Output:
%     NXi          - number of eigenfunctions with respect to SigmaXiP 
%                    retained for the analysis
%     NZeta        - number of eigenfunctions with respect to SigmaZetaP
%                    retained for the analysis
%     nu           - a (n*NXi+sumri*NZeta)x1 vector of scores
%     eveXi        - (pz*L0)xNXi Eigenfucntion 
%     eveZeta      - L0xNZeta Eigenfucntion 
%     evaXi        - NXi Eigenvalues 
%     evaZeta      - NZeta Eigenvalues 
%     SigmaPCA     - a p_z^2xL0 of optimal covariance function 
%     SigmaZetaPCA - L0x1 vector of variance matrix for Eta when s=t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude 

 L0=NoSetup(2);
 pz=NoSetup(5);

 [v1 d1]=eig(SigmaXiP);  
 dd1=diag(d1);   
 [dd1,ind1]=sort(dd1,'descend'); 
 [v2 d2]=eig(SigmaZetaP); 
 dd2=diag(d2);  
 [dd2,ind2]=sort(dd2,'descend');
 [~, d3]=eig(SigmaXiE);  
  dd3=diag(d3);  
 [dd3,~]=sort(dd3,'descend');
 
 dd11=dd1(dd1>0);
 dd22=dd2(dd2>0);
 dd33=dd3(dd3>0);
 tZeta=diag(SigmaZetaE);

 TotalVar=sum(dd11)+sum(dd22)+sum(dd33)+sum(tZeta(tZeta>0))/L0;

 if nargin == 6
      EV=0;
      NXi=0;
      NZeta=0;
      while(EV<ExpVar) 
          
     if (NXi+1<=length(dd11) && NZeta+1<=length(dd22))
          
          if (dd11(NXi+1)>= dd22(NZeta+1))
            NXi=NXi+1;  
          else
            NZeta=NZeta+1;   
          end     
     else   
         if (NXi+1>length(dd11))
            NZeta=NZeta+1;   
         end  
         if (NZeta+1>length(dd22))
            NXi=NXi+1;   
         end        
     end   
     
     EV=(sum(dd11(1:NXi))+sum(dd22(1:NZeta))+sum(dd33)+sum(tZeta(tZeta>0))/L0)/TotalVar;  
      end   
      
 end
 
 evaXi=dd1(1:NXi);
 evaZeta=dd2(1:NZeta);
 eveXi=v1(:,ind1(1:NXi)); 
 eveZeta=v2(:,ind2(1:NZeta));  

temp1=repmat(eveXi(:), 1,pz); 
temp11=temp1';
temp4=reshape(temp11,pz^2, L0*NXi); 

temp2=reshape(eveXi, pz,L0*NXi); 
temp3=repmat(temp2, pz,1);   

temp6=repmat(evaXi,1,L0)';
temp7=repmat(temp6(:)',pz^2,1); 
sigmaXiPCA=sum(reshape(temp7.*(temp4.*temp3), pz^2, L0, NXi),3);

temp8=repmat(evaZeta,1,L0)';  
temp9=eveZeta.*eveZeta;   
SigmaZetaPCA=sum(temp8.*temp9,2);   
  