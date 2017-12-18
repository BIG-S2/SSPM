function [SigmaZetaP,SigmaZetad,SigmaXiP,SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,ResYdesign)
% Input:
%     NoSetup      - col vector of [n, L0, p,sumri] 
%                    n = number of subjects 
%                    L0 = number of location of fiber tract 
%                    p = number of covariates 
%                    sumri = total number of measurements
%                    p_z = number of covariates correspoinding to random
%                    effects
%     Indicator    - nx1 cell, where  Indicator{i}=k1:k2 is a r_ix1
%                    vector, contains the index for each subject i in array
%                    Ydesign.r_i = number of time points for subject i.
%                    -----very important
%     arclength    - L0x1 col vector of the arclength from one end to the 
%                    other end
%     Zdesign      - (\sum_i=1^n ri) x p_z normalized design matrix
%     ResYdesign   - ((\sum_i=1^n ri)xL0)   matrix 
%     SigmaXi      - pzxpzxL0 covariance matrix for Xi
%     SigmaZeta    - L0xL0 covariance matrix for zeta
%%
% Output:
%     SigmaXiP     - (pzxL0) x (pzxL0) covariance matrix for Xi
%     SigmaZetaP   - L0xL0 covariance matrix for zeta
%     SigmaXid     - pzxL0 x xpzxL0 covariance matrix for Xi
%     SigmaZetad   - L0x1 covariance matrix for zeta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude 

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(5);
sumri=NoSetup(4);
sumII=0;
for nii=1:n
sumII=sumII+size(Indicator{nii},1)^2;
end
a2=1/sumri;

Index=(1:sumII)';
Indicator1=cell(n,1);
Indicator2=cell(n,1);
st=0;
for nii=1:n
    ri=size(Indicator{nii},1);
    st=st+1;
    ed=st-1+ri^2;
    Index1=Index(st:ed);
    Indicator2{nii}=Index1(1:(ri+1):ri^2); 
   
    Index1(1:(ri+1):ri^2)=[];
    Indicator1{nii}=Index1;        
    if (ri==1)
    Indicator1{nii}=[];    
    end
    st=ed;
end
Ind2=cell2mat(Indicator2);

p1=p*(p+1)/2;
M=zeros(p1,p^2);
k=0;
for i=1:p
    for j=1:p
        k=k+1;
        k1=0;
        for i1=1:p
            for j1=i1:p
                k1=k1+1;
                if (i1==i &&j1==j)||(i1==j &&j1==i)
                    M(k1,k)=1;
                end
                
            end
        end
    end
end


M2=zeros(p^2,p1);
k=0;
for i=1:p
    for j=1:p
        k=k+1;
        k1=0;
        for i1=1:p
            for j1=1:i1
                k1=k1+1;
                if ((i1==i&&j1==j)|| (i1==j && j1==i))
                    M2(k,k1)=1;
                end
                
            end
        end
    end
end

    Tempee=cell(1,n);  
    for nii=1:n      
        ri=size(Indicator{nii},1);
        tempRY=ResYdesign(Indicator{nii},:);   
        Tempee1=repmat(tempRY(:),1,ri)'; 
        Tempee11=reshape(Tempee1, ri^2, L0);     
        Tempee111=repmat(Tempee11, L0,1);  
        Tempee2=repmat(tempRY,ri,1);   
        Tempee22=repmat(Tempee2(:), 1,L0); 
       Tempee3=Tempee111.*Tempee22;    
       Tempee4=reshape(Tempee3,  ri^2, L0^2); 
       Tempee{nii}=Tempee4';   
    end
    TTT=cell2mat(Tempee);
    A2=reshape(sum(TTT(:,Ind2),2)*a2, L0,L0);    
   
TempXX1=cell(1,n);
TempXX2=cell(1,n);
TempXX1e=cell(1,n);
TempXX2e=cell(1,n);
for nii=1:n
     ri=size(Indicator{nii},1);
    TempX=Zdesign(Indicator{nii},:)';         
    
    TempX1=repmat(TempX(:),1,p)';            
    TempX11=reshape(TempX1, p^2,ri);        
    TempX111=repmat(TempX11,ri,1);        
    
    TempX2=repmat(TempX',1,p)';              
    TempX22=repmat(TempX2(:),1,ri);           
    
   TempXX=TempX111.*TempX22;               
   TTempXX1=reshape(TempXX, p^2, ri^2);
   TempXX1{1,nii}=TTempXX1;                                    
  
   MXX=M*TTempXX1;         
   TempXX1e{1,nii}=MXX;      
    
    TempXX1l=reshape(repmat(TempXX(:),1,p^2)', p^2*p^2,ri^2);  
    TempXX1r=repmat(reshape(TempXX, p^2, ri^2),p^2,1);     
    TempXX2{1,nii}=TempXX1l.*TempXX1r;      
    
    TempXX1le=reshape(repmat(MXX(:),1,p1)', p1*p1,ri^2);  
    TempXX1re=repmat(reshape(MXX, p1, ri^2),p1,1);      
    TempXX2e{1,nii}=TempXX1le.*TempXX1re;     
end
SSS=cell2mat(TempXX1);
Bn2=sum(SSS(:,Ind2),2);     
SSS1=cell2mat(TempXX2);
Dn= pinv(reshape(sum(SSS1,2), p^2,p^2));

SSSe=cell2mat(TempXX1e);
Bn2e=sum(SSSe(:,Ind2),2);    
SSS1e=cell2mat(TempXX2e);
Dne= pinv(reshape(sum(SSS1e,2), p1,p1));

tempA=SSS*TTT'; 
tS2=cell(L0,L0);
kk=0;
 for L0ii1=1:L0
        for L0ii2=1:L0
   kk=kk+1;   
         tS2{L0ii1,L0ii2}=tempA(:,kk);
        end
 end
A=cell2mat(tS2);


TTTe=TTT(1:(L0+1):L0^2,:);
tempAe=SSSe*TTTe';  
tS2e=cell(L0,1);
 for L0ii1=1:L0
         tS2e{L0ii1,1}=tempAe(:,L0ii1);
 end
Ae=cell2mat(tS2e); 

SigmaZetaP=((1-a2*(Bn2'*Dn*Bn2))^(-1))*(A2-kron(eye(L0), a2*Bn2'*Dn)*A); 
vecMSigmaXiP=kron(eye(L0), Dn)* (A-kron(SigmaZetaP,Bn2)); 
 
vecMSigmaXiP1=reshape(vecMSigmaXiP, p,p,L0, L0);  
vecMSigmaXiP1=permute(vecMSigmaXiP1, [2 1 3,4]);
vecMSigmaXiP2=permute(vecMSigmaXiP1, [1 3 2 4]);
SigmaXiP=reshape(vecMSigmaXiP2, p*L0, p*L0);

SigmaZetaPe=((1-a2*(Bn2e'*Dne*Bn2e))^(-1))*(diag(A2)-kron(eye(L0), a2*Bn2e'*Dne)*Ae); 
vecMSigmaXiPe=reshape(kron(eye(L0), Dne)* (Ae-kron(SigmaZetaPe,Bn2e)),p1,L0);

 for L0ii1=1:L0
 SigmaXiP(((L0ii1-1)*p+1):(L0ii1)*p,((L0ii1-1)*p+1):(L0ii1)*p)=reshape(M2*vecMSigmaXiPe(:,L0ii1),p,p);  
 SigmaZetaP(L0ii1,L0ii1)=SigmaZetaPe(L0ii1,1);
 end

SigmaZetad=diag(SigmaZetaP)';

SigmaXid=zeros(p^2,L0);
 for L0ii1=1:L0
 SigmaXid(:,L0ii1)=reshape(SigmaXiP(((L0ii1-1)*p+1):(L0ii1)*p,((L0ii1-1)*p+1):(L0ii1)*p),p^2,1);       
 end
