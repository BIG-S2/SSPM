
OutputImageMatrixFilesFullName=fullfile('C:\Documents and Settings\yuaih\My Documents\MATLAB\OUTPUT','OutputFilenames.txt'); 
fid=fopen(OutputImageMatrixFilesFullName);
dimX=256;
dimY=256;
dimZ=198;
Matr_temp=zeros(dimX,dimY,dimZ);
 tline = fgetl(fid);
 mm=1;
 while ischar(tline)
   fidinside=fopen(tline);
   if(fidinside~=-1)
   for kk=1:dimZ
       for ii=1:dimX
          Matr_temp(ii,:,kk)=fread(fidinside,dimY);
       end
   end
   fclose(fidinside);
   
   fprintf('%dth output images is being transformed...\n',mm);
   
   [Path,Name]=fileparts(tline);
   imName=fullfile(Path,Name);
   IM=make_nii(Matr_temp);
   save_nii(IM,imName);
   delete (tline);  
   tline = fgetl(fid);
   mm=mm+1;
   end
 end
 fclose(fid);    
 delete (OutputImageMatrixFilesFullName);