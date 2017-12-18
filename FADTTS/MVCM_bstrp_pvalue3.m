function [Gpval] = MVCM_bstrp_pvalue3(NoSetup,arclength,Xdesign,Ydesign,efitBetas1,InvSigmats,mh,Cdesign,B0vector,Gstat,GG,kstr) 
% 
% MVCM_bstrp_pvalue is to implement Zhu's (2010) method of calculating p value using bootstrap generated sample in MVCM
%
% Input:
%     NoSetup       - col vector of [n, L0, p, m] 
%                     n = sample size 
%                     L0 = number of location of fiber tract 
%                     p = number of covariates 
%                     m = number of features 
%     arclength     - col vector of the arclength from one end to the other  end
%     efitBetas1    - a 2*p x L0 x m matrix of estimated coefficients 
%     Xdesign       - a n x p normalized design matrix
%     Ydesign       - a n x L0 x m matrix of fiber bundle diffusion properties 
%     InvSigmats    - a 2p x 2p x L0 x m matrix
%     Cdesign       - a r x mp matrix for characterizing the r linear constriants among mp parameters 
%     B0vector      - a r x L0 vector for hypothesis testing
%     mh            - a m x 1 optimal bandwidth 
%     GG            - number of bootstrap
%     Gstat         - an over all test statitsics
%     kstr          - kernel function
% Output:
%     Gpval         - an over all p value
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%     [efitBetas,InvSigmats,efitYdesign] = MVCM_lpks_wb1(NoSetup,arclength,Xdesign,Ydesign,mh,kstr)
% before you use MVCM_bstrp_pvalue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% April 4, 2010 @ AA
%     

if nargin<14, 
  kstr='exp(-.5*t.^2)';
end  

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
m=NoSetup(4);

H0efitBetas1=efitBetas1*0.0;
AA=kron(eye(p),[1 0]);
BB=AA;
for mii=2:m
    BB=blkdiag(BB,AA);
end
CdesignNew=Cdesign*BB;

for L0ii=1:L0
    VecBeta10=squeeze(efitBetas1(:,L0ii,:));
    InvSigs=InvSigmats(:,:,L0ii,1);
    for mii=2:m
        InvSigs=blkdiag(InvSigs,InvSigmats(:,:,L0ii,mii));
    end
    VecBeta1=VecBeta10(:)-InvSigs*CdesignNew'*pinv(CdesignNew*InvSigs*CdesignNew')*(CdesignNew*VecBeta10(:)-B0vector(:,L0ii));
    H0efitBetas1(:,L0ii,:)=reshape(VecBeta1,[2*p,m]);
end

H0efitBetas=reshape(AA*H0efitBetas1(:,:),[p,L0,m]);
H0efitYdesign=reshape(Xdesign*H0efitBetas(:,:),[n,L0,m]);
H0ResYdesign=Ydesign-H0efitYdesign;
[H0ResEtas,H0efitEtas]=MVCM_sif2(arclength,H0ResYdesign);
        
        Gpval=0;
        Gstatvec=zeros(1,GG);
        Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';
        
        Tempmat=zeros(n,L0,2*p,L0,m);
        for mii=1:m
            h=mh(mii);
            Tmat=Tmat0/h;
            t=Tmat;
            Kmat=eval(kstr)/h;
            for nii=1:n
                for L0ii=1:L0
                    TempA=repmat(Xdesign(nii,:),[1,L0]);
                    TempB=kron(Tmat(:,L0ii)',ones(1,p));
                    TempC=[TempA;TempA.*TempB];
                    TempD=reshape(TempC(:),[2*p,L0]);
                    TempDD=permute(TempD,[2 1]);
                    TempEE=TempDD.*repmat(Kmat(:,L0ii),[1,2*p]);
                    Tempmat(nii,:,:,L0ii,mii)=TempEE;
                end
            end
        end

        for gii=1:GG
            
            [SimYdesign]=MVCM_grs(H0efitBetas,H0efitEtas,H0ResEtas,Xdesign);
                                 
            SimefitBetas=zeros(p,L0,m);
            SimResYdesign=SimYdesign*0.0;
            for mii=1:m
                AA=repmat(SimYdesign(:,:,mii),[1 1 2*p]);
                TempBeta=zeros(2*p,L0);
                for L0ii=1:L0
                    BB=Tempmat(:,:,:,L0ii,mii);
                    KXYZmat=sum(reshape(sum(AA(:,:).*BB(:,:)),[L0,2*p]))';
                    TempBeta(:,L0ii)=InvSigmats(:,:,L0ii,mii)*KXYZmat;
                end
                TempBeta1=kron(eye(p),[1 0])*TempBeta;
                SimefitBetas(:,:,mii)=TempBeta1;
                SimResYdesign(:,:,mii)=SimYdesign(:,:,mii)-Xdesign*TempBeta1;
            end

            [SimResEtas,SimefitEtas,SimeSigEta]=MVCM_sif2(arclength,SimResYdesign);
            
            [SimGstat]=MVCM_bstrp_stat(arclength,Xdesign,SimefitBetas,SimeSigEta,Cdesign,B0vector);
                        
            Gstatvec(1,gii)=SimGstat;
            Gpval=Gpval+(SimGstat>=Gstat);
        end
        Gpval=Gpval/GG;

               
end