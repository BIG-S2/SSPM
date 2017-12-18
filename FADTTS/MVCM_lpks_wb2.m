function [efitBetas3] = MVCM_lpks_wb2(NoSetup,arclength,Xdesign,Ydesign,mh,kstr)
% 
% MVCM_lpks_wb1 is to implement Zhu's (2010) method of local polynomial kernel smoothing (order = 3) with preselected bandwidth in MVCM
%
% Input:
%     NoSetup      - col vector of [n, L0, p, m] 
%                    n = sample size 
%                    L0 = number of location of fiber tract 
%                    p = number of covariates 
%                    m = number of features 
%     arclength    - col vector of the arclength from one end to the other end
%     Xdesign      - a n x p normalized design matrix. 
%     Ydesign      - a n x L0 x m matrix of fiber bundle diffusion properties 
%     mh           - a 1 x m vector of optimal bandwidth
%     kstr         - Kernel function  
% Output:
%     efitBetas3    - a 4*p x L0 x m matrix of estimated coefficients and dreivative of coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%    function [NoSetup, arclength, Xdesign,Ydesign]=MVCM_read(tractdata, designdata,  diffusionFiles, nofeatures, featuresname)  
% before you use MVCM_lpks_wb2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 27, 2010 @ AA
%     

if nargin<6, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
m=NoSetup(4);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';

efitBetas3=zeros(4*p,L0,m);

for mii=1:m
  
    h=mh(mii); 
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;

    Sigmat=zeros(4*p,4*p,L0);
    Tempmat0=zeros(4*p,4*p,L0);
    Tempmat3=zeros(L0,4*p,L0);
    Tempmat2=zeros(L0,4*p,L0);
    for nii=1:n
        for L0ii=1:L0
            TempA=repmat(Xdesign(nii,:),[1,L0]);
            TempB=kron(Tmat(:,L0ii)',ones(1,p));
            TempC=[TempA;TempA.*TempB;TempA.*TempB.^2;TempA.*TempB.^3];
            TempD=reshape(TempC(:),[4*p,L0]);
            TempDD=permute(TempD,[2 1]);
            TempEE=TempDD.*repmat(Kmat(:,L0ii),[1,4*p]);
            Tempmat2(:,:,L0ii)=TempEE.*repmat(Ydesign(nii,:,mii)',[1,4*p]);
            Tempmat0(:,:,L0ii)=TempEE'*TempDD;
        end
        Sigmat=Sigmat+Tempmat0;
        Tempmat3=Tempmat3+Tempmat2;
    end
    
    KXYZmat=reshape(sum(Tempmat3(:,:)),[4*p,L0]);
  
    for L0ii=1:L0
        efitBetas3(:,L0ii,mii)=pinv(Sigmat(:,:,L0ii))*KXYZmat(:,L0ii);
    end
   
end

end

