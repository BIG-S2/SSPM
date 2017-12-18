  % This is the main function for MAGEE method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Before running this program, please enter the required parameters,
%filenames and input image and output image directories bellow
 tic
DimSPD=1;   %Dimension of response
 num_Subjects=5;   %Number of subjects
 num_cov=4;   %Number of covariates for each dimension of response
 MaxTimePoints=5;    %Maximum number of time points among all subjects
 No_Row=1;   %Number of rows of constraint matrix
 num_iteration=6;   %It is suggested set to be between 5 and 20.
 
 ExampleDir='F:\SSPM-V2.1\MAGEE_example\';
 InputFileFolderDir=fullfile(ExampleDir,'InputFileFolder\');
 TimePointsFname='TimePointFilename.txt';
 DesignMatrixFname='DesignMatrixName.txt';
 MaskFname='Mask.hdr';
 ConsMatrixFname='linearConstraintMatrixName.txt';
 InputImageFn='InputImageFilename.txt';
 
 InputImageFilesDirectory=fullfile(ExampleDir,'\InputImages');
 OutputFilesDirectory=fullfile(ExampleDir,'\OutputImages');  % The output images directory.
 if exist(OutputFilesDirectory,'dir')~=7
    mkdir(OutputFilesDirectory);
 end       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %You do not need to chage anything bellow.  
 
 TimePointsFname=strcat(InputFileFolderDir,'\',TimePointsFname);
 DesignMatrixFname=strcat(InputFileFolderDir,'\',DesignMatrixFname);
 MaskFname=strcat(InputFileFolderDir,'\',MaskFname);
 ConsMatrixFname=strcat(InputFileFolderDir,'\',ConsMatrixFname);
 InputImageFn=strcat(InputFileFolderDir,'\',InputImageFn);
 
fprintf('Dimension of response:                    %d\n',DimSPD);
fprintf('Number of sujects:                        %d\n',num_Subjects);
fprintf('Number of covariates:                     %d\n',num_cov);
fprintf('Maximum number of time points:            %d\n',MaxTimePoints); 
fprintf('Number of rows of constraint matrix:      %d\n',No_Row);
fprintf('Number of iterations:                     %d\n\n',num_iteration);
fprintf('Filename containing time points:\n%s\n',TimePointsFname);
fprintf('Filename containing design matrix:\n%s\n',DesignMatrixFname);
fprintf('Mask image name: \n%s\n',MaskFname);
fprintf('Filename containing constraint matrix:\n%s\n',ConsMatrixFname);
fprintf('Filename containing input images:\n%s\n\n',InputImageFn);
fprintf('Input image file directory:\n%s\n',InputImageFilesDirectory);
fprintf('Output image file directory:\n%s\n\n',OutputFilesDirectory);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Mask_im=load_nii(MaskFname);
  Mask_matrix=Mask_im.img;
  [dimX,dimY,dimZ]=size(Mask_matrix);
  voxels=dimX*dimY*dimZ;
  [path,name,ext] = fileparts(MaskFname);  
  MaskImName=strcat(path,'\',name,'.img');
   
  fprintf('The mask image size : [%d, %d, %d]\n',dimX,dimY,dimZ);
  fprintf('Number of voxels: %d\n', voxels); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 InputImageFilesFullName = strcat(InputImageFilesDirectory,'\InputImageFilesname.txt');
 fprintf('InputImageFilesFullName:%s\n',InputImageFilesFullName);
 fid=fopen(InputImageFilesFullName,'wt');
 fid1=fopen(InputImageFn,'rt');
 tline = fgetl(fid1);
 ii=1;
 index=0;
while ischar(tline)
   str=strcat(InputImageFilesDirectory,'\',tline,'.img');
   fprintf(fid,'%s\n',str);
   im=load_nii(str);
   matr=im.img;
   [dimxr,dimyr,dimzr]=size(matr);
   fprintf('The %dth real image data size : [%d, %d, %d]\n',ii,dimxr,dimyr,dimzr);
   
   if dimX~=dimxr ||dimY~=dimyr || dimZ~=dimzr
      fprintf('The %dth real image data size does not match those of the mask image!\n',ii);
      index=1;
      break;
   end
   tline = fgetl(fid1);   
   ii=ii+1;   
end
fclose(fid1);
fclose(fid);
clear matr
if index==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% mex MAGEE.cpp  predef.cpp newdelete.cpp 
%%% compile mex file MAGEE main program
fprintf('\nData processing, please wait...\n');

  MAGEE(DimSPD,num_Subjects,num_cov,MaxTimePoints,TimePointsFname,DesignMatrixFname,...
      MaskImName,InputImageFilesFullName,No_Row,ConsMatrixFname,OutputFilesDirectory,dimX,dimY,dimZ,num_iteration);

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% transform output data into imge file
 fprintf('\nTransforming output data into imges,please wait...\n\n');
 OutputImageMatrixFilesFullName=strcat(OutputFilesDirectory,'\OutputFilenames.txt'); 
 fid=fopen(OutputImageMatrixFilesFullName,'rt');
 Matr_temp=zeros(dimX,dimY,dimZ,'single');
 tline = fgetl(fid);
 mm=1;
 while ischar(tline)
   fidinside=fopen(tline,'r');
   for kk=1:dimZ
       for jj=1:dimY
           Matr_temp(:,jj,kk)=fread(fidinside,dimX,'float32');
       end
   end
   fclose(fidinside);
   fprintf('%dth output images is being created...\n',mm); 
   [Path,Name]=fileparts(tline);
   imName=fullfile(Path,Name);  
   imName=strcat(imName,'.nii');        
   IM=make_nii(Matr_temp,[1,1,1],[1.0,1.0,1.0],'');
   save_nii(IM,imName);
   tline = fgetl(fid);
   mm=mm+1;
   Matr_temp=zeros(dimX,dimY,dimZ,'single');
 end
 fclose(fid); 
 
  str=strcat(OutputFilesDirectory,'\*.dat');
  delete(str)
  str=strcat(OutputFilesDirectory,'\OutputFilenames.txt');
  delete(str)
 

end
fprintf('++++++Your job is finished!!++++++++\n');

toc
