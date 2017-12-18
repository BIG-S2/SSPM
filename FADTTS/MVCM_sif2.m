function [ResEtas,efitEtas,eSigEta] = MVCM_sif2(arclength,ResYdesign,kstr) 
% 
% MVCM_sif is to implement Zhu's (2010) method of smoothing individual function without preselected bandwidth in MVCM
%
% Input:
%     arclength   - col vector of the arclength from one end to the other end
%     ResYdesign  - a n x L0 x m matrix of difference between fiber bundle diffusion properties and fitted fiber bundle diffusion properties
%     kstr        - Kernel function  
% Output:
%     ResEtas     - a n x L0 x m matrix of of difference between ResYdesign and fitted eta
%     efitEtas    - a n x L0 x m matrix of estimated eta
%     eSigEta     - a m x m x L0 x L0 matrix of covariance matrix of etas NOTE: only diagonal matrices are correct!!!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Please run  
%     [efitBetas,InvSigmats,efitYdesign] = MVCM_lpks_wb1(NoSetup,arclength,Xdesign,Ydesign,mh,kstr)
% before you use MVCM_sif2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% March 28, 2010 @ AA
%     

if nargin<3, 
  kstr='exp(-.5*t.^2)';
end 

[n L0 m]=size(ResYdesign);
xrange=max(arclength)-min(arclength);
nh=floor(max(30, L0/2));
hmin=1.0*xrange/L0; 
hmax=xrange/8; 
vh=logspace(log10(hmin),log10(hmax),nh);
Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';
GCVs=zeros(nh,m);

for nhii=1:nh
    h=vh(nhii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    S0=ones(L0,1)*sum(Kmat); 
    S1=ones(L0,1)*sum(Kmat.*Tmat0); 
    S2=ones(L0,1)*sum(Kmat.*Tmat0.^2);    
    Smat=Kmat.*(S2-S1.*Tmat0)./(S2.*S0-S1.^2);
    Smat=Smat';
    for mii=1:m
        GCVs(nhii,mii)=sum(sum(((ResYdesign(:,:,mii))'-Smat*(ResYdesign(:,:,mii))').^2))/(1-trace(Smat)/L0)^2/n;
    end
end

[MinGCV,flag]=min(GCVs);
mh=vh(flag);

efitEtas=zeros(n,L0,m);
%efitEta0=zeros(L0,m);

for mii=1:m
    h=mh(mii);
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    S0=ones(L0,1)*sum(Kmat); 
    S1=ones(L0,1)*sum(Kmat.*Tmat0); 
    S2=ones(L0,1)*sum(Kmat.*Tmat0.^2);    
    Smat=Kmat.*(S2-S1.*Tmat0)./(S2.*S0-S1.^2);
    efitEtas(:,:,mii)=ResYdesign(:,:,mii)*Smat;
    %efitEta0(:,mii)=squeeze(sum(efitEtas(:,:,mii)))/n;
end

ResEtas=ResYdesign-efitEtas;

%efitEta0All=efitEta0';
eSigEta=zeros(m,m,L0,L0);
for L0ii=1:L0
    %for L0jj=1:L0
        %eSigEta(:,:,L0ii,L0jj)=(squeeze(efitEtas(:,L0ii,:))-kron(ones(n,1),efitEta0All(:,L0ii)'))'*(squeeze(efitEtas(:,L0jj,:))-kron(ones(n,1),efitEta0All(:,L0jj)'))/(n-m);
        eSigEta(:,:,L0ii,L0ii)=(squeeze(efitEtas(:,L0ii,:)))'*squeeze(efitEtas(:,L0ii,:))/n;
        %eSigEta(:,:,L0ii,L0ii)=(squeeze(efitEtas(:,L0ii,:))-kron(ones(n,1),efitEta0All(:,L0ii)'))'*(squeeze(efitEtas(:,L0ii,:))-kron(ones(n,1),efitEta0All(:,L0ii)'))/(n-m);
    %end
end

end