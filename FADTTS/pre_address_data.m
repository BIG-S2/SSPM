
%%% Loading and Transferring data   %%% 

filepath='H:/PDF/Fibertract/Version_0812/Code/Data/Right_InternalCapsule_Cerebellar_Statistics';

%filepath2='H:/PDF/Fibertract/FRATS_Document_1015/GUI/FRACTS_GUI';
filepath2='H:\PDF\Fibertract\FADTTS_Document_Nov23\GUI';

%filepath2='D:/Linglong/Study/PDF/UNC/Fibertract/FRATS_Document_1015/GUI/FRACTS_GUI';

%filepath='D:/Linglong/Study/PDF/UNC/Fibertract/Version_0829/Code/Data/Splenium_Stats';

datapath='Right_InternalCapsule_Cerebellar_Statistics_LK2';

designdataname='Xdesigndata2.txt';
FAname='RightInternalCapsule-FA-LK2.txt';
MDname='RightInternalCapsule-MD-LK2.txt';
Lambda1name='RightInternalCapsule-Lambda1-LK2.txt';
Lambda2name='RightInternalCapsule-Lambda2-LK2.txt';
Lambda3name='RightInternalCapsule-Lambda3-LK2.txt';
tractname='tract.txt';

%%% Design data %%%
aa=sprintf('%s/%s',filepath,designdataname);
designdata1=load(aa);
designdata=designdata1(1:64,:);
aa=sprintf('%s/designdata.mat',filepath2);
save(aa,'designdata');

%%%  Combining FA, MD, lambda1, lambda2, and lambda3 together  %%%
%diffusionFiles=cell(5,1);
diffusionFiles=cell(2,1);
aa=sprintf('%s/%s/%s',filepath,datapath,FAname);
temp=load(aa);
diffusionFiles{1}=temp(:,2:65);
aa=sprintf('%s/%s/%s',filepath,datapath,MDname);
temp=load(aa);
diffusionFiles{2}=temp(:,2:65);
%aa=sprintf('%s/%s/%s',filepath,datapath,Lambda1name);
%temp=load(aa);
%diffusionFiles{3}=temp(:,2:65);
%aa=sprintf('%s/%s/%s',filepath,datapath,Lambda2name);
%temp=load(aa);
%diffusionFiles{4}=temp(:,2:65);
%aa=sprintf('%s/%s/%s',filepath,datapath,Lambda3name);
%temp=load(aa);
%diffusionFiles{5}=temp(:,2:65);
aa=sprintf('%s/diffusiondata.mat',filepath2);
save(aa,'diffusionFiles');

%%% tract data %%%
aa=sprintf('%s/%s',filepath,tractname);
tractdata=load(aa);
aa=sprintf('%s/tractdata.mat',filepath2);
save(aa,'tractdata');

%%% Cdesign matrix %%%
Cdesign=kron(eye(2),[0 1 0]);
aa=sprintf('%s/cdesign.mat',filepath2);
save(aa,'Cdesign');

%%% B0matrix %%%
B0matrix=zeros(2,75);
aa=sprintf('%s/b0vector.mat',filepath2);
save(aa,'B0matrix');

filepath2='H:/PDF/Fibertract/FRATS_Document_1015/GUI/FRACTS_GUI';
%%% Cdesign matrix %%%
Cdesign=[0 1 0];
aa=sprintf('%s/cdesign0.mat',filepath2);
save(aa,'Cdesign');

%%% B0matrix %%%
B0matrix=zeros(1,75);
aa=sprintf('%s/b0vector0.mat',filepath2);
save(aa,'B0matrix');
