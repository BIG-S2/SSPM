function [mh,GCVs,vh] = MVCM_lpks_wob(NoSetup,arclength,Xdesign,Ydesign,kstr)
% 
% MVCM_lpks_wob is to implement Zhu's (2010) method of local polynomial kernel smoothing (order = 1) without preselected bandwidth in MVCM
%
% Input:
%     NoSetup    - col vector of [n, L0, p, m] 
%                  n = sample size 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  m = number of features 
%     arclength  - col vector of the arclength from one end to the other end
%     Xdesign    - a n x p normalized design matrix. 
%     Ydesign    - a n x L0 x m matrix of fiber bundle diffusion properties 
%     kstr       - Kernel function  
% Output:
%     mh         - a 1 x m vector of optimal bandwidth
%     GCVs       - a nh x m matrix of GCVs
%     vh         - a 1 x nh vector of bandwidth
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%    function [NoSetup, arclength, Xdesign,Ydesign]=MVCM_read(tractdata, designdata,  diffusionFiles, nofeatures, featuresname)  
% before you use MVCM_lpks_wob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 27, 2010 @ AA
%     

if nargin<5, 
  kstr='exp(-.5*t.^2)';
end   

n=NoSetup(1);
L0=NoSetup(2);
p=NoSetup(3);
m=NoSetup(4);

xrange=max(arclength)-min(arclength);
nh=floor(max(30, L0/2));
hmin=1.0*xrange/L0; 
hmax=xrange/8; 
efitYdesign=Ydesign*0;
GCVs=zeros(nh,m);
vh=logspace(log10(hmin),log10(hmax),nh);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';


for nhii=1:nh
    h=vh(nhii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;

    Tempmat0=zeros(2*p,2*p,L0,n);
    Tempmat=zeros(n,L0,2*p,L0);
    for nii=1:n
        for L0ii=1:L0
            TempA=repmat(Xdesign(nii,:),[1,L0]);
            TempB=kron(Tmat(:,L0ii)',ones(1,p));
            TempC=[TempA;TempA.*TempB];
            TempD=reshape(TempC(:),[2*p,L0]);
            TempDD=permute(TempD,[2 1]);
            TempEE=TempDD.*repmat(Kmat(:,L0ii),[1,2*p]);
            Tempmat(nii,:,:,L0ii)=TempEE;
            Tempmat0(:,:,L0ii,nii)=TempEE'*TempDD;
        end
    end

    Tempmat1=zeros(n,2*p,L0,m);
    for mii=1:m        
        AA=repmat(squeeze(Ydesign(:,:,mii)),[1 1 2*p L0]);
        Tempmat10=Tempmat(:,:).*AA(:,:); 
        Tempmat11=reshape(Tempmat10,[n,L0,2*p,L0]);        
        Tempmat1(:,:,:,mii)=squeeze(sum(Tempmat11,2));        
    end

    for mii=1:m
        for L0ii=1:L0
            Tempmat2=sum(squeeze(Tempmat0(:,:,L0ii,:)),3);
            Tempmat3=sum(squeeze(Tempmat1(:,:,L0ii,mii)),1);
            for nii=1:n
                Tempmat22=Tempmat2-Tempmat0(:,:,L0ii,nii);
                Tempmat33=Tempmat3-squeeze(Tempmat1(nii,:,L0ii,mii));
                AA=pinv(Tempmat22)*Tempmat33';
                efitYdesign(nii,L0ii,mii)=Xdesign(nii,:)*kron(eye(p),[1 0])*AA;            
            end
        end
        GCVs(nhii,mii)=sum(sum((Ydesign(:,:,mii)-efitYdesign(:,:,mii)).^2))/(n*L0); 
    end

end

[MinGCV,flag]=min(GCVs);
mh=vh(flag);

end

