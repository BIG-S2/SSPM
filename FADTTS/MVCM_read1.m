function [Xdesign, Ydesign,scalediffusion]=MVCM_read1(designdata,respdata)  

% MVCM_read1 is to standarddize the Xdesign and Ydesign matrix
% Inputs:
%     designData: the text file containing covariates of interest. Please always include the intercept in the first column. 
%                 is a n x p matrix. 
%                 n denotes the number of subjects. 
%                 p denotes the number of covariates. 
%     respData:  n x L0 x m matrix of response.
% Output:
%     Xdesign    - a n x p normalized design matrix. 
%     Ydesign    - a n x L0 x m matrix  
%     scalediffusion: a m x 1 vector of scales for each properties 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
%    Copyright (c) H. T. Zhu  2009 
%    Modifiedy by L. L. Kong 2010

%% Normalize the design matrix     
 
Xdesign=designdata; 
if size(Xdesign, 1)< size(Xdesign, 2) 
   Xdesign=Xdesign'; 
end
for i=1:size(Xdesign, 2)  
   if std(Xdesign(:,i))>0 
      Xdesign(:,i)=(Xdesign(:,i)-mean(Xdesign(:,i)))/std(Xdesign(:,i));  
   end  
end 

%% Standardized all diffusion properties to a comparable scale 
Ydesign0=permute(respdata,[3 2 1]);
nofeatures=size(Ydesign0,1);

meanYdata=zeros([nofeatures, 1]); 
for i=1:nofeatures 
   meanYdata(i)=mean(mean(abs(Ydesign0(i,:,:)))); 
end 
maxMeanYdata=max(meanYdata); 
for i=1:nofeatures 
   Ydesign0(i,:,:)=Ydesign0(i,:,:)*floor(maxMeanYdata/meanYdata(i)); 
end 
scalediffusion=floor(meanYdata./maxMeanYdata); 

Ydesign=permute(Ydesign0,[3 2 1]);

end 

