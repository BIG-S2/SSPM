function [mh] = BiSmooth_LL_bs(NoSetup,arclength,Sigma,p,In,kstr)

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
%     In           - 0 0r 1 
%     kstr         - Kernel function 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mh           - bandwidth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude 
if nargin<6, 
  kstr='exp(-.5*t.^2)';
end   

L0=NoSetup(2);
xrange=max(arclength)-min(arclength);
nh=floor(L0/6);
hmin=1.0*xrange/L0; 

hmax=xrange/8; 
GCVs=zeros(nh,1);
vh=logspace(log10(hmin),log10(hmax),nh);

Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

SigmaK1=reshape(permute(reshape(Sigma,p,L0,p,L0), [1,3,2,4]),p^2,L0^2)';
 
 if In==0
 SigmaK1(1:(L0+1):L0^2,:)=0; 
 end
 
 SigmaK2=reshape(SigmaK1,L0, L0,p^2);
      
for nhii=1:nh
   
    h=vh(nhii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    Sigma1=cell(L0,L0);
    
    for L0ii1=1:L0
         for L0ii2=1:L0ii1  
             
          SigmaK22=SigmaK2;
          SigmaK22(L0ii1,L0ii2,p^2)=0;
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
  
        att=reshape(permute(reshape(repmat(TempT(:)',3,1),3,3,L0^2), [1 3 2]),3*L0^2,3);
          TempmatL1=TempKK1.*(repmat(TempT(:),1,3).*att); 
          TempmatL2=reshape(TempmatL1',3^2,L0^2);
          
          TempmatL22=reshape(TempmatL2, 3^2,L0,L0);
    
          TL=pinv(reshape(sum(TempmatL2,2)-squeeze(TempmatL22(:,L0ii1,L0ii2)),3,3));
          TR= (TL*TempmatR1)*reshape(SigmaK22,L0^2,p^2); 
          
          Sigma1{L0ii1,L0ii2}=reshape(TR(1,:),p,p);
          Sigma1{L0ii2,L0ii1}=reshape(TR(1,:),p,p)';
          end
    end
    
Sigma11=cell2mat(Sigma1);    
v=Sigma-Sigma11;
GCVs(nhii)=sqrt(v(:)'*v(:));   
end

[~,flag]=min(GCVs);
mh=vh(flag);
