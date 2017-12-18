function [Ydesign,scalediffusion]=MVCM_read2(respdata)  

% MVCM_read1 is to standarddize the Xdesign and Ydesign matrix
% Inputs:
%     respData:  n x L0 x m matrix of response.
% Output:
%     Ydesign    - a n x L0 x m matrix  
%     scalediffusion: a m x 1 vector of scales for each properties 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
%    Copyright (c) H. T. Zhu  2009 
%    Modifiedy by L. L. Kong 2010

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

