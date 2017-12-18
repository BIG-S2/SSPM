function [NoSetup, Xdesign,Zdesign]=FMPM_ReadNormalize(Xdesign,Zdesign,Indicator, arclength)  

% FMPM_ReadNormalize is to read the arclength, design matrix for functional
% mixed effects model and to normalize the design matrix.
% Inputs:
%      Xdesign     - (\sum_i=1^n ri) x p design matrix.
%      Zdesign     - (\sum_i=1^n ri) x p_z design matrix. 
%      Indicator   - nx1 cell, where  each cell is a r_ix1
%                    vector, contains the indices for each subject i
%      arclength   - L0x1 col vector of the arclength from one end to the 
%                    other end  
%
% Output:
%     NoSetup    - col vector of [n, L0, p, sumri]
%                  n = number of subjects 
%                  L0 = number of location of fiber tract 
%                  p = number of covariates 
%                  sumri = number of measurements for all n subjects  
%                  p_z=number of covariates corresponding to random effect
%     Xdesign    - a (\sum_i=1^n ri) x p normalized design matrix.
%     Zdesign    - a (\sum_i=1^n ri) x p_z normalized design matrix.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%    Copyright  by Y. Yuan 2014 at St. Jude

%  Set parameters and defaults according to number of input arguments
%
if nargin < 4;    
   exit(1); 
end;

NoSetup=zeros(4, 1);
NoSetup(1)=size(Indicator,1);   %n;
NoSetup(2)=size(arclength,1);   %L0;
NoSetup(3)=size(Xdesign,2);     %p;  including intercept;
NoSetup(4)=size(Xdesign,1);     %sumri;
NoSetup(5)=size(Zdesign,2);     %p_z

%% Normalize the design matrix   
for i=1:size(Xdesign, 2)  
   if std(Xdesign(:,i))>0 
      Xdesign(:,i)=(Xdesign(:,i)-mean(Xdesign(:,i)))/std(Xdesign(:,i));  
   end  
end 

%% Normalize the design matrix   
for i=1:size(Zdesign, 2)  
   if std(Zdesign(:,i))>0 
      Zdesign(:,i)=(Zdesign(:,i)-mean(Zdesign(:,i)))/std(Zdesign(:,i));  
   end  
end 
