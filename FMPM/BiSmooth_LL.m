 function [Sigma,SigmaE, SigmaE1] = BiSmooth_LL(NoSetup,arclength,Sigma, Sigmad1,p,mh,In,kstr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     NoSetup      - col vector of [n, L0, p,sumri] 
%                    n = number of subjects 
%                    L0 = number of location of fiber tract 
%                    p_x = number of covariates 
%                    sumri = total number of measurements
%                    p_z = number of covariates correspoinding to random
%                    effects
%     arclength    - L0x1 col vector of the arclength from one end to the 
%                    other end
%     Sigma        - (pxL0)x(pxL0)
%     Sigmad1      -  p^2xL0
%     In           - 0 0r 1 
%     kstr         - Kernel function 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sigma    -(pxL0)x(pxL0)
%   SigmaE   - pxp;
%   SigmaE1  -(pxL0)x(pxL0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude 

if nargin<8, 
  kstr='exp(-.5*t.^2)';
end   

Sigma0=Sigma;

L0=NoSetup(2);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';
    h=mh;
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    
 SigmaK1=reshape(permute(reshape(Sigma,p,L0,p,L0), [1,3,2,4]),p^2,L0^2)'; %L0^2*p^2
   
    TempmatR=zeros(3,L0^2,L0,L0);
    TempmatL=zeros(3^2,L0^2,L0,L0);
    Sigma1=cell(L0,L0);
    Sigmad=zeros(p^2,L0);
    for L0ii1=1:L0
         for L0ii2=1:L0ii1            
           TempK1=repmat(Kmat(:,L0ii1),L0,1);
           TempK2=repmat(Kmat(:,L0ii2),1,L0)';   
           TempK=TempK1.*TempK2(:); 
            if In==0
           TempK(1:(L0+1):L0^2)=0;
            end
           
           TempKK=repmat(TempK', 3,1); 
           TempKK1=repmat(reshape(TempKK,3*L0^2,1),1,3);  

           TempT1=repmat(Tmat(:,L0ii1),L0,1);
           TempT11=TempT1';   
           TempT2=repmat(Tmat(:,L0ii2),1,L0)';
           TempT22=TempT2(:)';
           TempT=[ones(1,L0^2);TempT11; TempT22]; 
           TempmatR1=TempKK.*TempT;       
           TempmatR(:,:,L0ii1,L0ii2)=TempmatR1;
           
          att=reshape(permute(reshape(repmat(TempT(:)',3,1),3,3,L0^2), [1 3 2]),3*L0^2,3);
          TempmatL1=TempKK1.*(repmat(TempT(:),1,3).*att); 
          TempmatL2=reshape(TempmatL1',3^2,L0^2);
          
          TempmatL(:,:,L0ii1,L0ii2)=TempmatL2;
          
          TL=pinv(reshape(sum(TempmatL2,2),3,3));
          TR= (TL*TempmatR1)*SigmaK1; 
          Sigma1{L0ii1,L0ii2}=reshape(TR(1,:),p,p);
          Sigma1{L0ii2,L0ii1}=reshape(TR(1,:),p,p)';
          if (L0ii2==L0ii1)
          Sigmad(:,L0ii1)=TR(1,:)'; 
          end
          
          end
    end
    
    Sigma=cell2mat(Sigma1);
    SigmaE=reshape(sum(Sigmad1-Sigmad,2)/L0,p,p);
    SigmaE1=Sigma0-Sigma;
  