function [ResEtas,efitEtas] = FMPM_LL_If_wob(NoSetup,arclength,ResYdesign,kstr) 
% 
% FMPM_LL_If_wob is to implement method of smoothing individual function 
% without preselected bandwidth in FMPM
%
% Input:
%     NoSetup    - col vector of [n, L0, p, sumri,p_z] 
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = total number of measurements
%                  p_z = number of covariates correspoinding to random
%                  effects
%     arclength   - col vector of the arclength from one end to the other end
%     ResYdesign  - a sumri x L0 matrix of difference between fiber bundle
%                   diffusion properties and fitted fiber bundle diffusion 
%                   properties
%     kstr        - Kernel function  
% Output:
%     ResEtas     - a sumri x L0 matrix of of difference between ResYdesign
%                   and fitted eta
%     efitEtas    - a sumri x L0 matrix of estimated eta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright  by Y. Yuan 2014 at St. Jude    

if nargin<4, 
  kstr='exp(-.5*t.^2)';
end 

L0=NoSetup(2);

xrange=max(arclength)-min(arclength);
nh=floor(L0/2);
hmin=1.0*xrange/L0; 
hmax=xrange/8; 

GCVs=zeros(nh,1);
vh=logspace(log10(hmin),log10(hmax),nh);

Tmat0=arclength*ones(1,L0)-ones(L0,1)*arclength';


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
  
   GCVs(nhii)=sum(sum(((ResYdesign)'-Smat*(ResYdesign)').^2))/(1-trace(Smat)/L0)^2;
  
end

[~,flag]=min(GCVs);
mh=vh(flag);

    h=mh;
    Tmat=Tmat0/h;
    t=Tmat;
    Kmat=eval(kstr)/h;
    S0=ones(L0,1)*sum(Kmat); 
    S1=ones(L0,1)*sum(Kmat.*Tmat0); 
    S2=ones(L0,1)*sum(Kmat.*Tmat0.^2);    
    Smat=Kmat.*(S2-S1.*Tmat0)./(S2.*S0-S1.^2);
    efitEtas=ResYdesign*Smat;
   
ResEtas=ResYdesign-efitEtas;
