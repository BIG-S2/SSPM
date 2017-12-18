function [NoSetup, arclength, Xdesign, Ydesign, scalediffusion]=MVCM_read(tractdata, designdata,  diffusionFiles, nofeatures)  

% MVCM_read is to read the fiber coordinate data, design matrix for functional linear model, and fiber bundle diffusion properties 
%   Can use first  3, 4 or 5 arguments.
% Inputs:
%     tractData: the text file containing (x, y, z) coordinates of all locations on a given fiber tract. 
%                The data set should start from one end to the other end. 
%                 tractData is a L0 x 3 matrix
%                 L0 denotes the number of locations.  
%                 3 denotes the three coordinates.  
%     designData: the text file containing covariates of interest. Please always include the intercept in the first column. 
%                 is a n x p matrix. 
%                 n denotes the number of subjects. 
%                 p denotes the number of covariates. 
%     diffusionFiles: 
%                 diffusionFiles is a mx1 cell containing the names of all fiber diffusion properties files.  
%                 Each fiber bundle diffusion properties should contain a L_0 x n matrix. 
%                      Rows correspond to the rows in tractData, while columns correspond to the columns in designData. 
%     nofeatures:  the number of diffusion properties that we want to joint model. 
%     scalediffusion: a m x 1 vector of scales for each properties 
% Output:
%     NoSetup    - col vector of [n, L0, p, m] 
%     arclength  -  col vector of the arclength from one end to the other end.
%     Xdesign    - a n x p normalized design matrix. 
%     Ydesign    - a n x L0 x m matrix  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
%    Copyright (c) H. T. Zhu  2009 
%    Modifiedy by L. L. Kong 2010



%  Set parameters and defaults according to number of input arguments
%
if nargin < 3 || nargin > 5;    %  less than 3 argument inputs,  quit 
   exit(1); 
end;

if nargin == 3 ;    %  3 arguments input,  m=1 
  nofeatures = 1;
end ;

%% Initinalizing matrices

NoSetup=zeros(4, 1); 


%% Calculating arc length 

tractCON=tractdata;  
if size(tractCON, 2)>3 
   tractCON=tractCON';
end 

NoSetup(2)=size(tractCON, 1);  %L0 
NoSetup(4)=nofeatures;  %L0 
arclength=zeros(size(tractCON,1), 1); 
temp1=tractCON-tractCON;
temp1(1,:)=tractCON(1,:);
temp1(2:size(tractCON,1),:)=tractCON(1:(size(tractCON,1)-1),:); 
temp1=temp1-tractCON; 
for i=1:size(tractCON, 1) 
    arclength(i)=norm(temp1(i,:)); 
    if i>1
       arclength(i)=arclength(i)+arclength(i-1); 
    end 
end 



%% Normalize the design matrix     
 
Xdesign=designdata; 
if size(Xdesign, 1)< size(Xdesign, 2) 
   Xdesign=Xdesign'; 
end
NoSetup(1)=size(Xdesign, 1);  % n 
NoSetup(3)=size(Xdesign, 2);  % p  
for i=1:size(Xdesign, 2)  
   if std(Xdesign(:,i))>0 
      Xdesign(:,i)=(Xdesign(:,i)-mean(Xdesign(:,i)))/std(Xdesign(:,i));  
   end  
end 

%% Combine all diffusion properties from the select fiber tract 

Ydesign0=zeros([nofeatures, NoSetup(2), NoSetup(1)]); 

for i=1:nofeatures 
 tempt=diffusionFiles{i}; 
 if NoSetup(2)==size(tempt, 1) 
  Ydesign0(i,:,:)=tempt(:,1:NoSetup(1));
 else 
   if NoSetup(2)==size(tempt, 2) 
     tempt=tempt';
     Ydesign0(i,:,:)=tempt(:,1:NoSetup(1));
   else 
     exit(1); 
   end 
 end 
end 

%% Standardized all diffusion properties to a comparable scale 

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

