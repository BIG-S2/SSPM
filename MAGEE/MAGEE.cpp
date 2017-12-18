// MARM.cpp : Defines the entry point for the console application.
//
# include <cstdio> 
# include <cstdlib>  
# include <cmath>    
# include <fstream>
# include <string.h>
#include "predef.h"
#include "newdelete.h"
#include "mex.h"

int  DimSPD, DimXX, DimYY, htNoM,  htNoR, Ncov, NcovDim, NcovVar;
int  NOmaxTIMEs,  NopointSur, NOtraG, Nsample, TotalImg, NumRep;

unsigned int congrval, tausval;

/*****************External Functions************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
	int i, j, k, i1, ii, jj, kk,  dii, djj, dkk, dmm, dll,  dimX, dimY, dimZ;
	int num_subjects, num_cov, noRow;
	int seed, SEEDY[12];
	int pointINDX;
 	double hMARM; 
	FILE  *dat2, *dat3, *runlog;


	//****************************************************
	//*    Get Output Files Directory from Input File
	//****************************************************

	char *OutputDir;
	int buflen;
	buflen = mxGetN(prhs[10])+1;
    OutputDir = (char *) mxCalloc(buflen,sizeof(char));       
    mxGetString(prhs[10], OutputDir , buflen); 
	// OutputDir:Output Files Directory  
	
    char* logname=new char[200];
	dll=sprintf(logname, "%s/MAGEE_run_log.txt",OutputDir); 
	runlog=fopen(logname, "w");
	delete[] logname;

	char* logMesg=new char[500];
    dll=sprintf(logMesg, "Output directory:%s",OutputDir); 
	fprintf(runlog,"%s\n",logMesg);

	//***********************************
	// Initinalizing random numbers 
	//***********************************
	
	seed=17;
    setSEED(SEEDY, seed);
    init1(SEEDY);

	//***********************************
	// Enter the setup parameters  
	//***********************************
    DimSPD=(int)mxGetScalar(prhs[0]);  //Dimension of response;
	dll=sprintf(logMesg,"Dimension of response that mexfunction read:%d\n",DimSPD);
	fprintf(runlog,"%s",logMesg);
	

	DimYY=DimSPD;
	num_subjects=(int)mxGetScalar(prhs[1]); // number of subjects;
    Nsample=num_subjects;
	
	dll=sprintf(logMesg, "number of subjects that mexfunction read:%d",num_subjects); 
	fprintf(runlog,"%s\n",logMesg);
		
	num_cov=(int)mxGetScalar(prhs[2]);   //number of covairiates;
   	dll=sprintf(logMesg, "number of covairiates that mexfunction read:%d",num_cov);
	fprintf(runlog,"%s\n",logMesg);
		
	NOmaxTIMEs=(int)mxGetScalar(prhs[3]);   //maximum number of timepoint for each suject;
    dll=sprintf(logMesg, "maximum number of timepoint for each suject that mexfunction read:%d",NOmaxTIMEs);
	fprintf(runlog,"%s\n",logMesg);

	hMARM = 1.15;
	  
    int*      allDimTime=ivector(1, num_subjects+1);  
    int*      CCID=ivector(0, DimYY*num_cov); 
	float**   chi2ref=matrix(1, 100, 1, 9);  
	float**   ExactTime=matrix(1, num_subjects+1, 1, NOmaxTIMEs+1); 
    int*      indxMi=ivector(0, num_subjects); // reading the number of time points insides  HT 
    int*      indxMi2=ivector(0, num_subjects); // reading the number of time points insides  HT 
	float**   INRR=matrix(0, DimYY*num_cov, 0, DimYY*num_cov);
    float**   INRRnew=matrix(0, DimYY*num_cov, 0, DimYY*num_cov);
	double**  INRRnewD=dmatrix(1, DimYY*num_cov, 1, DimYY*num_cov);
	int*      Irank=ivector(0, DimYY*num_cov);
	int*      mi=ivector(0,num_subjects+1);
    float**   midRR=matrix(0, DimYY*num_cov, 0, DimYY*num_cov);
	float**   RR=matrix(0, DimYY*num_cov-1, 0, DimYY*num_cov-1); // linear constraints: Wald test
    float*    rr00=vector(0, DimYY*num_cov);
    double**  TXX=dmatrix(1, num_cov*DimYY, 1, num_cov*DimYY);	
	
	float**  TXXtemp1=matrix(1, num_cov*DimYY, 1, num_cov*DimYY);
    float**  TXXtemp2=matrix(1, num_cov*DimYY, 1, num_cov*DimYY);


   //****************************************************  
   // Enter chi-square quantiles                       // 
   //****************************************************/

     chisqQuantile(chi2ref); 

   //****************************************************  
   // Enter file name containing Time points            // 
   //****************************************************

	double minTime=0.0, maxTime=0.0; 
    
	
    char *tempo, *tempo_1, *tempo_2;  
   
	buflen = mxGetN(prhs[4]) + 1;
    tempo_1 = (char *) mxCalloc(buflen,sizeof(char));//file name containing Time points
     
    mxGetString(prhs[4], tempo_1, buflen);  
	dll=sprintf(logMesg, "file name containing Time points that mexfunction read:%s",tempo_1);
	fprintf(runlog,"%s\n",logMesg);
   	
	 if((dat3=fopen(tempo_1, "r"))==NULL){
		mexPrintf("Cannot open the file name containing Time points.\n");
		dll=sprintf(logMesg, "Cannot open the file name containing Time points.");
	    fprintf(runlog,"%s\n",logMesg);
   	   
		return;
	 }  
    mxFree(tempo_1);

	TotalImg=0; 
	maxTime=0.0; 
	fprintf(runlog, "\n");
	fprintf(runlog, "allDimTime:\n");
	for (dii = 1; dii <= num_subjects; dii++){
		fscanf(dat3, "%d", &allDimTime[dii]);
		dll = sprintf(logMesg, "\n %d ", allDimTime[dii]);
		fprintf(runlog, "%s\n", logMesg);
		mi[dii] = (int)(allDimTime[dii]);
		TotalImg += allDimTime[dii];
		if (mi[dii]>maxTime)
			maxTime = mi[dii];
	}
	
	NOmaxTIMEs=maxTime; 
	maxTime=0; 

   mexPrintf("\n");

    for(dii=1; dii<=Nsample; dii++){
		mexPrintf("\n %d ", allDimTime[dii]);
		fprintf(runlog, "\n");
		fprintf(runlog, "allDimTime:\n");
		dll = sprintf(logMesg, "\n %d ", allDimTime[dii]);
		fprintf(runlog, "%s\n", logMesg);
		for (djj = 1; djj <= allDimTime[dii]; djj++){
			fscanf(dat3, "%f", &ExactTime[dii][djj]);
			dll = sprintf(logMesg, "%f ", ExactTime[dii][djj]);
			fprintf(runlog, "%s\n", logMesg);
            mexPrintf("%f ", ExactTime[dii][djj]);  
		  if(minTime>ExactTime[dii][djj])
			  minTime=ExactTime[dii][djj];
		  if(maxTime<ExactTime[dii][djj])
			   maxTime=ExactTime[dii][djj]; 
	   }
   	} 
     
	 mexPrintf("\n\n");
	 dll = sprintf(logMesg, "minTime:%f; maxTime:%f  ", minTime, maxTime);
	 fprintf(runlog, "%s\n", logMesg);
	 fprintf(runlog, "\n\n");
	 printf("\n");
	 fclose(dat3);

	// Standardize the time points // 

   for(dii=1; dii<=Nsample; dii++) 
	  for(djj=1; djj<=allDimTime[dii]; djj++)
		 if(maxTime-minTime>0)  
	        ExactTime[dii][djj]=(ExactTime[dii][djj]-minTime)*4.0/(maxTime-minTime*1.0);
		 else
            ExactTime[dii][djj]=0.0; 

	int allDimMax=0;
	for(dii=1; dii<=Nsample; dii++){
		if (allDimTime[dii]>allDimMax)
			allDimMax=allDimTime[dii];
	}

	indxMi[0]=0; 
	indxMi2[0]=0; 
	indxMi[1]=allDimTime[1]; 
    indxMi2[1]=(int)(allDimTime[1]*(allDimTime[1]-1)*0.5); 
    mexPrintf("%d %d %d\n", 1, indxMi[1], indxMi2[1]);
	for(dii=2; dii<=Nsample; dii++){
		indxMi[dii]=allDimTime[dii]; 
		indxMi2[dii]=(int)(allDimTime[dii]*(allDimTime[dii]-1)*0.5);  
        indxMi[dii]+=indxMi[dii-1]; 
		indxMi2[dii]+=indxMi2[dii-1];
		mexPrintf("%d %d %d\n", dii, indxMi[dii], indxMi2[dii]); 
	} 

     mexPrintf("\n"); 
	//******************************  
   // Enter design matrix 
   //******************************

	float**   x0x= matrix(0, TotalImg-1, 0, num_cov-1);
	
    buflen = mxGetN(prhs[5])+1;
    tempo_2 = (char *) mxCalloc(buflen,sizeof(char)); //design matrix file name     
    mxGetString(prhs[5], tempo_2, buflen); 

	dll=sprintf(logMesg, "design matrix file name:%s",tempo_2); 
	fprintf(runlog,"%s\n",logMesg);
	
   	
	 if((dat2=fopen(tempo_2, "r"))==NULL){
		mexPrintf("Cannot open the design matrix file name.\n");
		dll=sprintf(logMesg, "Cannot open the design matrix file name.");
	    fprintf(runlog,"%s\n",logMesg);
		return;//////////
	 }
	 for (dii = 0; dii<TotalImg; dii++)
		 for (djj = 0; djj<num_cov; djj++){
			 fscanf(dat2, "%f", &x0x[dii][djj]);
		 }
	 fclose(dat2);
    mxFree(tempo_2);

	 
	 //******************************  
	// Enter mask image file  
	//******************************
	
	 int NopointSurNew=0;
	 char *maskTempo;

	buflen = mxGetN(prhs[6])+1;
    maskTempo = (char *) mxCalloc(buflen,sizeof(char)); 
    mxGetString(prhs[6], maskTempo, buflen); 
    //maskTempo: mask image name
	dll=sprintf(logMesg, "mask image name:%s",maskTempo); 
	fprintf(runlog,"%s\n",logMesg);

    dimX=(int)mxGetScalar(prhs[11]); 
	dimY=(int)mxGetScalar(prhs[12]); 
	dimZ=(int)mxGetScalar(prhs[13]); 

	int voxels=dimX*dimY*dimZ;
	float*** mask_tmp = fmatrix3(0, dimZ-1, 0, dimX-1, 0, dimY-1);
	readAnalyzeBody_char((char*)maskTempo, dimZ, dimX, dimY, mask_tmp);
	mxFree(maskTempo);
   // NopointSur = 0;
    int index=1;
    for(k=0; k<dimZ; k++)
	    for (j = 0; j<dimY; j++)
    	   for(i=0; i<dimX; i++)
    		{
				if(mask_tmp[k][i][j] > 0.1)
					NopointSurNew++;
				index++;
			}
   	//****************************************************
	//*    Get Input Image Matrix Filename from InputImageFile
	//****************************************************

    buflen = mxGetN(prhs[7])+1;
    tempo = (char *) mxCalloc(buflen,sizeof(char));       
    mxGetString(prhs[7], tempo , buflen);  
	//tempo: Input Image Matrix Filename
	dll=sprintf(logMesg, "Input Image Matrix Filename:%s",tempo); 
	fprintf(runlog,"%s\n",logMesg);
   
	//***************************
	//*    Input Image  Data
	
	int**    cord=imatrix(1, NopointSurNew, 1, 3);
	int*     flag=ivector(1,dimX*dimY*dimZ);
	float**  ImageData=matrix(1, NopointSurNew, 1, TotalImg*DimYY); 
	float*** tensors = fmatrix3(0, dimZ-1, 0, dimX-1, 0, dimY-1);
	double*  Yresponse=dvector(1, TotalImg*DimYY);
	double     maxYvalue=0.0; 	
   	
    FILE* fp;
	int len;	
	
	if((fp = fopen(tempo, "r"))==NULL){
		mexPrintf("Cannot open the input image matrix file\n");
		mexPrintf("input image matrix file name: %s\n",tempo);
		dll=sprintf(logMesg, "input image matrix file name: %s",tempo);
	    fprintf(runlog,"%s\n",logMesg);
		mexCallMATLAB(0,NULL,0,NULL,"pause");//////
		exit(1);//////////
	  }
	mxFree(tempo);
	dmm=0;
	mexPrintf("\nTotalImg=%d, DimSPD=%d\n",TotalImg,DimSPD);
	mexPrintf("\n");
	dll=sprintf(logMesg, "\nTotalImg=%d, DimSPD=%d\n",TotalImg,DimSPD);
	fprintf(runlog,"%s\n",logMesg);

    char* eachTensorFile=(char *) mxCalloc(200,sizeof(char)); 
	
	for(i=0; i<TotalImg; i++){
  	   for(djj=0;djj<DimSPD;djj++){	
		   mexPrintf(" Reading (%d %d)th image\n", i, djj); 
		   fgets(eachTensorFile,200,fp);           
			len=strlen(eachTensorFile);
			char* eachTensorFile_new=(char *) mxCalloc(200,sizeof(char)); 
			strncpy(eachTensorFile_new,eachTensorFile,len-1);
			readAnalyzeDTI(eachTensorFile_new, dimZ, dimX, dimY, 1, tensors); 
			mxFree(eachTensorFile_new);
			index=1;
			for(k=0; k<dimZ; k++)
			   for(j=0;j<dimY;j++)
    			  for(dii=0; dii<dimX; dii++){
					   if(mask_tmp[k][dii][j]> 0.1){  				
						ImageData[index][i*DimSPD+djj+1]=tensors[k][dii][j];
						if(ImageData[index][i*DimSPD+djj+1]>maxYvalue)
						maxYvalue=ImageData[index][i*DimSPD+djj+1];   			       
						index++;
					}
				}  
			dll = sprintf(logMesg, "maxYvalue of image data = %1f\n", maxYvalue);
			fprintf(runlog, "%s", logMesg);
         }//djj
	}// i 


	mxFree(eachTensorFile);	
	free_fmatrix3(tensors, 0, 1*dimZ-1,  0, dimX-1, 0, dimY-1); 
	fclose(fp);
	
	mexPrintf("\n"); 
//*****************************
	//*    Input Constraint Matrix
	//*****************************
    
	noRow=(int)mxGetScalar(prhs[8]);//    noRow=2;  
    dll=sprintf(logMesg, "number of rows of constraint matrix that mexfunction read:%d",noRow); 
	fprintf(runlog,"%s\n",logMesg);	
	//number of rows of constraint matrix
	
	for(dii=0; dii<num_cov*DimYY; dii++)
		for(djj=0; djj<num_cov*DimYY; djj++){
  		    INRR[dii][djj]=0.0;
			midRR[dii][djj]=0.0;
		 }
	
	if(noRow>DimYY*num_cov-1)
		noRow=DimYY*num_cov-1;
	for(ii=0; ii<=DimYY*num_cov-1; ii++){
		for(jj=0; jj<=DimYY*num_cov-1; jj++)
			RR[ii][jj]=0.0;
        rr00[ii]=0.0;
	}
	
	char* fileNameG;	
    
    buflen = mxGetN(prhs[9])+1;
    fileNameG = (char *) mxCalloc(buflen,sizeof(char));       
    mxGetString(prhs[9], fileNameG , buflen); 
	//fileNameG:constraint matrix name

	dll=sprintf(logMesg, "constraint matrix name:%s",fileNameG); 
	fprintf(runlog,"%s\n",logMesg);
	
	if((dat2=fopen(fileNameG,"r"))==NULL){
       mexPrintf("Cannot open the save file\n");
	   mexPrintf("Input Constraint Matrix flie name: %s.\n",fileNameG);
	   dll=sprintf(logMesg, "Input Constraint Matrix flie name: %s.",fileNameG);
	   fprintf(runlog,"%s\n",logMesg);
		//mexCallMATLAB(0,NULL,0,NULL,"pause");
		exit(1);//////////
    }
	
	for(dii=0; dii<noRow; dii++){
		for(djj=0; djj<num_cov*DimYY; djj++){
		  fscanf(dat2, "%f", &RR[dii][djj]);
	      mexPrintf(" %f", RR[dii][djj]);
		}
		mexPrintf("\n");
	}
	mexPrintf("\n");
	for(dii=0; dii<noRow; dii++){
        fscanf(dat2, "%f", &rr00[dii]);
		mexPrintf("%f  ",rr00[dii]);
	}
	mexPrintf("\n");
	for(dii=0; dii<noRow; dii++){
		fscanf(dat2, "%d", &CCID[dii]);
		mexPrintf("%d ",CCID[dii]);
		CCID[dii]=CCID[dii]-1;
		midRR[0][CCID[dii]]=1;
	}
	mexPrintf("\n\n");
	fclose(dat2);

	for(dii=0, djj=noRow; dii<num_cov*DimYY; dii++)
		 if(midRR[0][dii]<1.0){
            CCID[djj]=dii;
            djj++;
		 }
    for(dii=0; dii<noRow; dii++)
		for(djj=0; djj<noRow; djj++)
		   INRR[dii][djj]=RR[dii][CCID[djj]];
    for(dii=0; dii<noRow; dii++){
		for(djj=0; djj<noRow; djj++)
			mexPrintf("%f",INRR[dii][djj]);
		mexPrintf("\n");
	 }
	mexPrintf("\n");

	ivYS(INRR, 0, noRow-1);
	for(dii=0; dii<noRow; dii++)
		for(djj=0; djj<num_cov*DimYY; djj++){
			midRR[dii][djj]=0.0;
			for(dkk=0; dkk<noRow; dkk++)
                midRR[dii][djj]-=INRR[dii][dkk]*RR[dkk][djj];
		 }
     for(dii=0; dii<noRow; dii++)
		 for(djj=0; djj<noRow; djj++)
             midRR[dii][CCID[djj]]=INRR[dii][djj];
     for(dii=noRow; dii<num_cov*DimYY; dii++)
		  midRR[dii][CCID[dii]]=1.0;

      //this generates midRR=[R_1^{-1}, -R_1^{-1}R_2]
      for(dii=0; dii<num_cov*DimYY; dii++)
		 for(djj=0; djj<num_cov*DimYY; djj++)
			 INRR[dii][djj]=midRR[dii][djj];

	  for(dii=0; dii<num_cov*DimYY; dii++)
		  for(djj=0; djj<num_cov*DimYY; djj++)
			  INRRnew[dii][djj]=0.0;
	  for(dii=0; dii<noRow; dii++)
          for(djj=0; djj<num_cov*DimYY; djj++)
			 INRRnew[dii][djj]=RR[dii][djj];
	  for(dii=noRow; dii<num_cov*DimYY; dii++)
		  INRRnew[dii][CCID[dii]]=1.0;
	  for(dii=1; dii<=num_cov*DimYY; dii++) 
		  for(djj=1; djj<=num_cov*DimYY; djj++) 
		     INRRnewD[dii][djj]=INRRnew[dii-1][djj-1];

     double  ht=DnewIVrank(INRRnewD, Irank, num_cov*DimYY);   
	/* for(dii=1; dii<=num_cov*DimYY; dii++){
		  for(djj=1; djj<=num_cov*DimYY; djj++) 
		       mexPrintf("%lf ", INRRnewD[dii][djj]);  
		  mexPrintf("\n"); 
	  } */
	 mexPrintf("\n");
	 
	
	int p=num_cov;
	int q=1;//variance

	Ncov=num_cov*DimYY;
	NcovVar=q;
	Nsample=num_subjects;

    mxFree(fileNameG);

	free_ivector(CCID, 0, DimYY*num_cov); 
	free_matrix(INRR, 0, DimYY*num_cov, 0, DimYY*num_cov);
	free_matrix(INRRnew, 0, DimYY*num_cov, 0, DimYY*num_cov);
    free_dmatrix(INRRnewD, 1, DimYY*num_cov, 1, DimYY*num_cov);
	free_ivector(Irank, 0, DimYY*num_cov);
	free_matrix(midRR, 0, DimYY*num_cov, 0, DimYY*num_cov);   

	

	//*****************************
	//*   Creat Design Matrix
	//*****************************	
	
	int totalobs, obsforonemax, mimax;
	double**   designXX = dmatrix(1, TotalImg*DimYY, 1, Ncov);
	
 
    for(dii=1; dii<=TotalImg*DimYY; dii++)
		for(dkk=1; dkk<=Ncov; dkk++) 
            designXX[dii][dkk]=0.0; 

    for(dii=1; dii<=TotalImg; dii++)
		for(dkk=1; dkk<=DimYY; dkk++) 
		   for(djj=1; djj<=num_cov; djj++) 
			designXX[(dii-1)*DimSPD+dkk][(dkk-1)*num_cov+djj]=(float)(x0x[dii-1][djj-1]*1.0);

    for(dii=1; dii<=Ncov; dii++)
	   for(djj=1; djj<=Ncov; djj++)
		   TXX[dii][djj]=0.0; 

	for(dii=1; dii<=TotalImg*DimYY; dii++)
		for(djj=1; djj<=Ncov; djj++) 
			for(dkk=1; dkk<=Ncov; dkk++)
			 TXX[djj][dkk]+=designXX[dii][djj]*designXX[dii][dkk];

    for(dii=1; dii<=Ncov; dii++)
		for(djj=1; djj<=Ncov; djj++)
			TXXtemp1[dii][djj]=(float)TXX[dii][djj]; 
 
		
	ABSmatrixNhalf(TXXtemp1, TXXtemp2, Ncov, -1.0); 
	
	for(dii=1; dii<=Ncov; dii++)
		for(djj=1; djj<=Ncov; djj++)
			TXX[dii][djj]=TXXtemp2[dii][djj]; 

	matrixprint(TXX, Ncov, Ncov); 

	
	totalobs=0;
	mimax=0;

	for(dii=1; dii<=num_subjects; dii++){
		if(mimax<=allDimTime[dii])
		   mimax=allDimTime[dii];
	  }
    totalobs=TotalImg*DimYY;
   	obsforonemax=mimax*DimYY;

	free_matrix(TXXtemp1, 1, num_cov*DimYY, 1, num_cov*DimYY);
	free_matrix(TXXtemp2, 1, num_cov*DimYY, 1, num_cov*DimYY);

	//**********************************************************
	//*  First step of MAGEE: Carry out voxel-wise Wald test.  
	//**********************************************************
	int newindex;
	int Nrowtotal;
    double tempValue, tempMean; 
	double**  allsurface=dmatrix(1,NopointSurNew,1,TotalImg*DimYY);
	double*   beta=dvector(1,Ncov+2*DimYY);
	float**   BrainTheta=matrix(1, NopointSurNew, 1, Ncov+2*DimYY);    
	float***  BrainThetaVar=fmatrix3(1,NopointSurNew,1,Ncov,1,Ncov);
	int*      flag2=ivector(1,voxels);
	double**  linearbeta=dmatrix(1,NopointSurNew,1,Ncov);
	int*      Nrow=ivector(1,Nsample);
	double*   residualli=dvector(1, totalobs);
	double**  varGEE=dmatrix(1, Ncov, 1,Ncov);
	float**** varMatrix=fmatrix4(1, NopointSurNew,1,Nsample,1, obsforonemax,1, obsforonemax);
	//float**   VarRes=matrix(1,NopointSurNew,1,1);
	double*   waldvaluetemp=dvector(1,NopointSurNew);
	mimax=0;
	Nrowtotal=0;
	for (dii=1;dii<=Nsample;dii++){
		Nrow[dii]=mi[dii]*DimYY;
		if (mi[dii]>mimax)
			mimax=mi[dii];
 		Nrowtotal+=Nrow[dii];
	}
	/*for (dii = 1; dii <= NopointSurNew; dii++)
		VarRes[dii][1] = 0.0;*/
  
	newindex=1; 
    pointINDX=1;
	index=0;
	mexPrintf("max=%lf\n", maxYvalue);
	
	for(k=0; k<dimZ; k++)
	     for (j = 0; j<dimY; j++)
		    for(dii=0; dii<dimX; dii++){
				index++;
				flag[index]=0;
				if(mask_tmp[k][dii][j]> 0.1){
					flag[index]=1; 
				//	if(index%1000==0)
				//			mexPrintf("%d\n", index); 
	        	     		tempValue=0.0; tempMean=0.0;  
	        		for(i=0; i<TotalImg; i++)
	        		{
 						for(ii=1; ii<=DimYY; ii++){ 
                            Yresponse[i*DimYY+ii]=ImageData[pointINDX][i*DimYY+ii]*10.0/maxYvalue;
                            tempValue+=Yresponse[i*DimYY+ii]*Yresponse[i*DimYY+ii]; 
							tempMean+=Yresponse[i*DimYY+ii]; 
						}	                   
					}
                    tempValue=tempValue/(TotalImg*1.0); 
                    tempMean=tempMean/(TotalImg*1.0); 
                    tempValue=sqrt(tempValue-tempMean*tempMean);  
				    if(fabs(tempValue)<=0.00001){
							flag[index]=0;
						//	mexPrintf("%lf ", tempValue); 
						}
	
					
					if(flag[index]>=1)
					{
						flag2[k*dimX*dimY + j*dimX + dii + 1] = newindex;

					 for(i=0; i<TotalImg; i++)
						  for(ii=1; ii<=DimYY; ii++){ 
						    allsurface[newindex][i*DimYY+ii]=ImageData[pointINDX][i*DimYY+ii];
						  }
					  cord[newindex][1]=k+1;
					  cord[newindex][2]=dii+1;
					  cord[newindex][3]=j+1;
		              GEEestimatesHongtu(beta,   residualli, designXX, TXX, allsurface[newindex], DimYY, allDimTime, Ncov, q, varMatrix[newindex], ExactTime,   indxMi,   indxMi2);
					  CalVarGEEHongtu(designXX, residualli, DimSPD, allDimTime, Ncov, beta,  varMatrix[newindex],  varGEE);
						
					  for (dkk = 1; dkk <= Ncov + 2 * DimYY; dkk++)
					      if (newindex == 1000) {
						  dll = sprintf(logMesg, "beta[%4d]=%lf\n", dkk, beta[dkk]);
						  fprintf(runlog, "%s", logMesg);
						  }
					  
					  for(dkk=1; dkk<=Ncov+2*DimYY; dkk++){
 					      BrainTheta[newindex][dkk]=beta[dkk];
						  if (newindex == 1000) {
							  dll = sprintf(logMesg, "BrainTheta[1000][%4d]=%lf\n", dkk, beta[dkk]);
							  fprintf(runlog, "%s", logMesg);
						  }
						}
				      for(dkk=1;dkk<=Ncov;dkk++)
					     for(dll=1;dll<=Ncov;dll++)
					  	    BrainThetaVar[newindex][dkk][dll]=varGEE[dkk][dll];
                      waldvaluetemp[newindex]=fabs(waldtest(BrainTheta[newindex],varGEE,RR,rr00,noRow));
					  if (waldvaluetemp[newindex] > 800.0){
						  for (dkk = 1; dkk <= Ncov + 2 * DimYY; dkk++)
							  mexPrintf(" %lf\n", beta[dkk]);
						  
						  for (dkk = 1; dkk <= Ncov; dkk++){
							  for (dll = 1; dll <= Ncov; dll++)
								  mexPrintf(" %lf ", varGEE[dkk][dll]);
							  mexPrintf("\n");
						  }
					  }
					  if (newindex % 1000 == 0){
						  mexPrintf("Wald[%4d]=%lf\n", newindex, waldvaluetemp[newindex]);
						  dll = sprintf(logMesg, "Wald[%4d]=%lf\n", newindex, waldvaluetemp[newindex]);
						  fprintf(runlog, "%s", logMesg);
					  }
					  newindex++;
					}
					pointINDX++;
				}//if mask

			} //j
			            
			free_fmatrix3(mask_tmp, 0, dimZ-1, 0, dimX-1, 0, dimY-1);
			free_dmatrix(TXX,1, num_cov*DimYY, 1, num_cov*DimYY);
			free_ivector(allDimTime,1, num_subjects+1); 
            free_matrix(ExactTime, 1, num_subjects+1, 1, NOmaxTIMEs+1);
            free_ivector(indxMi, 0, num_subjects); 
			free_ivector(indxMi2, 0, num_subjects); 	

            mexPrintf("\n END OF STEP 1: Carry out voxel-wise Wald test.\n");
			dll=sprintf(logMesg, "END OF STEP 1: Carry out voxel-wise Wald test.");
	        fprintf(runlog,"%s\n",logMesg);
  //******************************************************************************
  //   MARM  Step2:   Define the neighborhoods at each voxel
  //******************************************************************************

    NopointSur=newindex-1; //From here NopointSur change to newindex-1 

    int dllmax;
	double  hmax, Distance;
	int     neighindex, hnum;
	int     tempzl,tempzh,tempxl,tempxh,tempyl,tempyh,dnn;
	int*    noneighs=ivector(1,NopointSur);
	int**   neighID=imatrix(1,NopointSur,1,300);
	float** distance=matrix(1,NopointSur,1,300);	

	hnum=(int)mxGetScalar(prhs[14]); // number of iterations
	
	mexPrintf("\n pointINDX=%d, NopointSurNew=%d, NopointSur=%d, newindex=%d\n", pointINDX, NopointSurNew, NopointSur, newindex);
	
	dll=sprintf(logMesg, "pointINDX=%d, NopointSurNew=%d, NopointSur=%d, newindex=%d", pointINDX, NopointSurNew, NopointSur, newindex);
	fprintf(runlog,"%s\n",logMesg);

	hmax=pow(hMARM, hnum*1.0);


	index=1;
	pointINDX=1;

    for(dii=1;dii<=NopointSur;dii++)
		for(djj=1;djj<=Ncov;djj++)
		    linearbeta[dii][djj]=BrainTheta[dii][djj];
	
	dllmax=0;
    for(kk=1; kk<=NopointSur; kk++)
 		    {				
            tempzl=dmax(1,cord[kk][1]-10);
			tempzh=dmin(dimZ,cord[kk][1]+10);
			tempxl=dmax(1,cord[kk][2]-10);
			tempxh=dmin(dimX,cord[kk][2]+10);
			tempyl=dmax(1,cord[kk][3]-10);
			tempyh=dmin(dimY,cord[kk][3]+10);
            dll=0;
			for (djj = tempzl; djj <= tempzh; djj++)
			  for (dmm = tempyl; dmm <= tempyh; dmm++)
			     for (dkk = tempxl; dkk <= tempxh; dkk++){
				   if (flag[(djj - 1)*dimX*dimY + (dmm - 1)*dimX + dkk] >= 1){
						Distance = 0.0;
						Distance = 1.0*(cord[kk][1] - djj)*(cord[kk][1] - djj) + 1.0*(cord[kk][2] - dkk)*(cord[kk][2] - dkk) + 1.0*(cord[kk][3] - dmm)*(cord[kk][3] - dmm);
						if (Distance <= hmax){
							dll++;
							distance[kk][dll] = Distance;
							neighindex = (djj - 1)*dimX*dimY + (dmm - 1)*dimX + dkk;
							neighID[kk][dll] = flag2[neighindex];
  	                     }//if
			     	}//if
		     }//dmm
		     noneighs[kk]=dll;
			 if(dll>dllmax)
				dllmax=dll;
			}
 
	free_dvector(Yresponse, 1, TotalImg*DimYY);
	free_ivector(flag2,1,voxels);
	free_imatrix(cord,1,NopointSurNew,1,3);
	free_matrix(ImageData,1, NopointSurNew, 1, TotalImg*DimYY); 

	mexPrintf("\n END OF STEP 2: Define the neighborhoods at each voxel\n");
	dll=sprintf(logMesg, "END OF STEP 2: Define the neighborhoods at each voxel");
	fprintf(runlog,"%s\n",logMesg);
	// *************************************************************
	//         Step 3: Adaptive and Weighted Testing and Estimation               
	// *************************************************************
	
	int  doo;
	double h,c0,u1,kloc,u2,kst,tempvalue,qvalue, stoppingRULE;
	float temphongtu;
	float*** BrainThetaVarOld=fmatrix3(1,NopointSur,1,Ncov,1,Ncov);
	double** chivalue=dmatrix(1,1,1,1);
	float**  correctedwaldpvalue=matrix(1,hnum,1,NopointSur);
	int*     countbig=ivector(1,hnum);
	int*     countnumber=ivector(1,NopointSur);
    double** diff=dmatrix(1,Ncov,1,1);
	double** diff3=dmatrix(1,Ncov,1,1);
	int*     done=ivector(1,NopointSur);
	double** h3beta=dmatrix(1,NopointSur,1,Ncov);
	float**  holdbeta=matrix(1,NopointSur,1,Ncov);
	float*** INVvarGEE3=fmatrix3(1,NopointSur,1,Ncov,1,Ncov);
	int** labelid=imatrix(1,NopointSur, 1,dllmax+1);
	float**  sortwaldpvalue=matrix(1,hnum,1,NopointSur);
	float** varGEE3D=matrix(1,Ncov,1,Ncov);
    float**  waldpvalue=matrix(1,hnum,1,NopointSur);
	double*  waldpvalueold=dvector(1,NopointSur);
	float**  waldvalue=matrix(1,hnum,1,NopointSur);
	double*  waldvalueold=dvector(1,NopointSur);
	double** weight=dmatrix(1,NopointSur,1,dllmax+1);
	

	stoppingRULE=chi2ref[Ncov][5]; 
	FILE  *outf, *datSResidual1t;
	
	h=0.0;
	qvalue=0.1;
	temphongtu=0.0;
	
    for(djj=1; djj<=NopointSur; djj++)
		done[djj]=0;
	for(dii=1;dii<=hnum;dii++)
		countbig[dii]=0;

	c0 = log(Nsample*1.0)*chi2ref[Ncov][1];
	int countnotemp;
	int countwrong = 0;

	char* outputfilename = new char[200];
	dll = sprintf(outputfilename, "%s/OutputFilenames.txt", OutputDir);
	outf = fopen(outputfilename, "w");
	delete[] outputfilename;
	
	char *outdataname = new char[200];

	for (dii=1;dii<= hnum;dii++){ 
		for(djj=1; djj<=NopointSur; djj++){
			doo=1;
			if (h==0.0) {
				weight[djj][doo]=1.0;
				labelid[djj][doo]=djj;
				countnumber[djj]=doo;
				CalMeanParaWeightGEECorrectedForVar(designXX, BrainTheta[djj], mi, DimSPD, allsurface, weight[djj],labelid[djj],countnumber[djj],varMatrix);
			}
			else{
				if(done[djj]==0){
					countnotemp=0;
					for(dkk=1;dkk<=noneighs[djj];dkk++)
						if(distance[djj][dkk]<=h) 
							countnotemp+=1;
					for(dkk=1;dkk<=noneighs[djj];dkk++){
						if(distance[djj][dkk]<=h){
							u1=distance[djj][dkk]/h;
							if (u1<=1)
								kloc=1-u1;
							else kloc = 0.0;
							for (dll = 1; dll <= Ncov; dll++)
								diff[dll][1] = holdbeta[djj][dll] - holdbeta[neighID[djj][dkk]][dll];

							u2=0.0;
 						       for(dmm=1;dmm<=Ncov;dmm++) 
 								 for(dnn=1;dnn<=Ncov;dnn++) 
									u2+=diff[dmm][1]*INVvarGEE3[djj][dmm][dnn]*diff[dnn][1];

						/*	if(u2<0.0){
								for(dmm=1; dmm<=Ncov; dmm++) 
									mexPrintf(" %lf ", diff[dmm][1]); 
								mexPrintf("\n voxel=[%d]  %lf \n",djj,   u2);  	
							}*/
							kst=exp(-fabs(u2)/c0);
							weight[djj][doo]=kst*kloc*1.0;
							labelid[djj][doo]=neighID[djj][dkk];
							doo++;
						}//if
					}//dkk: calculate weight			
					countnumber[djj]=doo-1;//Yimei May 18

					CalMeanParaWeightGEECorrectedForVar(designXX, BrainTheta[djj], mi, DimSPD, allsurface, weight[djj], labelid[djj], countnumber[djj], varMatrix);
					if(dii==3){
						for(dkk=1;dkk<=Ncov;dkk++)
							h3beta[djj][dkk]=BrainTheta[djj][dkk];
					}
					if(dii>=4){
						for(dll=1;dll<=Ncov;dll++){						
							diff3[dll][1]=BrainTheta[djj][dll]-h3beta[djj][dll];
						}
						chivalue[1][1]=0.0; 
						for(dll=1; dll<=Ncov; dll++)
							for(dkk=1; dkk<=Ncov; dkk++) 
								chivalue[1][1]+=diff3[dll][1]*INVvarGEE3[djj][dll][dkk]*diff3[dkk][1]; 
			//			mexPrintf("MARM(%d)=%lf\n", dii, chivalue[1][1]);
						if(chivalue[1][1]<0.0){
							for(dll=1; dll<=Ncov; dll++) 
								mexPrintf("%lf %lf", INVvarGEE3[djj][dll][dll], diff3[dll][1]); 
							mexPrintf("\n"); 	
						}
						if(fabs(chivalue[1][1])>stoppingRULE/sqrt(dii*1.0)) {   //chisq(2,0.8)
							done[djj]=1;
							for(dkk=1;dkk<=Ncov;dkk++)
								BrainTheta[djj][dkk]=holdbeta[djj][dkk];
						}//if chivalue
					}//if dii
				}// if done[djj]
				else {
					for(dkk=1;dkk<=Ncov;dkk++)
						BrainTheta[djj][dkk]=holdbeta[djj][dkk];
				}
			}//else	 if h^=0.0	        
		}//djj

		if (dii==1||dii==5||dii==hnum){       
			for(dmm=1; dmm<=Ncov; dmm++){  
				dll=sprintf(outdataname, "%s/BetaMA_h%d_P%d.dat",OutputDir, dii, dmm);
				fprintf(outf, "%s\n", outdataname);
				datSResidual1t=fopen(outdataname, "wb+");
				for(dkk=1, djj=1; dkk<=voxels; dkk++){
					if(flag[dkk]>=1)
					{
						fwrite(&BrainTheta[djj][dmm], sizeof(float), 1, datSResidual1t);
						djj+=1;
					}else
						fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
				}/* dkk */ 
				fclose(datSResidual1t);
                dll=sprintf(outdataname, "%s/BetaSTD_h%d_P%d.dat", OutputDir,dii, dmm);
				fprintf(outf, "%s\n", outdataname);
				datSResidual1t=fopen(outdataname, "wb+");
				for(dkk=1, djj=1; dkk<=voxels; dkk++){
					if(flag[dkk]>=1)
					{
						fwrite(&BrainThetaVar[djj][dmm][dmm], sizeof(float), 1, datSResidual1t);
						djj+=1;
					}else
						fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
				}/* dkk */ 
				fclose(datSResidual1t);
							 
			}/* dmm */ 
		}/* if */ 

		for (djj = 1; djj <= NopointSur; djj++){
			if (done[djj] == 0){
				CalVarWeightedGEECorrectedForVar(djj, allsurface, designXX, DimSPD, mi, BrainTheta, q, varGEE, weight[djj], labelid[djj], countnumber[djj], varMatrix);
				for (dll = 1; dll <= Ncov; dll++)
					for (dmm = 1; dmm <= Ncov; dmm++)
						BrainThetaVar[djj][dll][dmm] = varGEE[dll][dmm];
				waldvalue[dii][djj] = waldtest(BrainTheta[djj], varGEE, RR, rr00, noRow);
				tempvalue = fabs(waldvalue[dii][djj])*(Nsample - noRow) / (noRow*(Nsample - 1)*1.0);
				if (dii == 1 || dii == 5 || dii == hnum){
					waldpvalue[dii][djj] = betai((Nsample - noRow)*1.0 / 2.0, noRow*1.0 / 2.0, (Nsample - noRow) / ((Nsample - noRow)*1.0 + noRow*tempvalue*1.0));
					if (waldpvalue[dii][djj] <= 0.0){
						waldpvalue[dii][djj] = 0.0;
						//printf("Large=%f \n", waldvalue[dii][djj]); 
					}
					sortwaldpvalue[dii][djj] = waldpvalue[dii][djj];
					waldvalueold[djj] = waldvalue[dii][djj];
					waldpvalueold[djj] = waldpvalue[dii][djj];
				}
				else{
					waldvalue[dii][djj] = waldvalueold[djj];
					waldpvalue[dii][djj] = waldpvalueold[djj];
				}
			}// if done[djj]
			else {
				waldvalue[dii][djj] = waldvalueold[djj];
				waldpvalue[dii][djj] = waldpvalueold[djj];
				sortwaldpvalue[dii][djj] = waldpvalueold[djj];
				for (dll = 1; dll <= Ncov; dll++)
					for (dmm = 1; dmm <= Ncov; dmm++){
						BrainThetaVar[djj][dll][dmm] = BrainThetaVarOld[djj][dll][dmm];
						varGEE[dll][dmm] = BrainThetaVarOld[djj][dll][dmm];
				}
			}

		}//djj
		
		if(dii==1){
			for (djj=1;djj<=NopointSur;djj++){
			  for (dll=1;dll<=Ncov;dll++)
			 	 for (dmm=1;dmm<=Ncov;dmm++)
					 varGEE3D[dll][dmm]=(float) BrainThetaVar[djj][dll][dmm];	
			  ABSmatrixNhalf(varGEE3D, INVvarGEE3[djj], Ncov, -1.0); 
			}
		}//if dii==1
		
		for(djj=1;djj<=NopointSur;djj++){
			for(dkk=1;dkk<=Ncov;dkk++)
				holdbeta[djj][dkk]=BrainTheta[djj][dkk];
			for (dkk=1;dkk<=Ncov;dkk++)
				for(dll=1;dll<=Ncov;dll++)
					BrainThetaVarOld[djj][dkk][dll]=BrainThetaVar[djj][dkk][dll];
			correctedwaldpvalue[dii][djj]=waldpvalue[dii][djj]*NopointSur;
			if(correctedwaldpvalue[dii][djj]>=1.0)
				correctedwaldpvalue[dii][djj]=1.0;
		}//djj
		
		for (djj=1;djj<=NopointSur;djj++){
			if (fabs(waldpvalue[dii][djj])<=0.000) 
				waldpvalue[dii][djj]=100.0;
			else
				waldpvalue[dii][djj]=-log10(waldpvalue[dii][djj]*1.0);
			if (fabs(correctedwaldpvalue[dii][djj])<=0.0) 
				correctedwaldpvalue[dii][djj]=60.0;
			else
				correctedwaldpvalue[dii][djj]=-log10(correctedwaldpvalue[dii][djj]*1.0);
			if (correctedwaldpvalue[dii][djj]>2.0) countbig[dii]+=1;
				}//djj            
  	    if (h==0.0)
			h=1.15;
		else 
			h=h*1.15;
	}//dii

    for (dii=1;dii<=hnum;dii++)
		mexPrintf("\n countbig[dii]=%d",countbig[dii]);

	free_matrix(chi2ref,1, 100, 1, 9);
    free_ivector(mi, 0,num_subjects+1);
	free_matrix(RR, 0, DimYY*num_cov-1, 0, DimYY*num_cov-1); 

	mexPrintf("\n\n END OF STEP 3: Adaptive and Weighted Testing and Estimation\n");
	dll=sprintf(logMesg, "END OF STEP 3: Adaptive and Weighted Testing and Estimation");
	fprintf(runlog,"%s\n",logMesg);

// float temphongtu;
//**************************************
//*        Output Results             
//**************************************
	temphongtu=0.0;

	if (hnum <= 5) {
	
		for(dmm=1; dmm<=2; dmm++){     
			if(dmm==1)
				i1=1;
			if(dmm==2)
				i1=hnum; 		
			dll=sprintf(outdataname, "%s/BwaldMA_h%d.dat",OutputDir, i1);
			fprintf(outf, "%s\n", outdataname);
			datSResidual1t=fopen(outdataname, "wb+");
			for(dkk=1, djj=1; dkk<=voxels; dkk++){
				if(flag[dkk]>=1)
				{
					fwrite(&waldvalue[i1][djj], sizeof(float), 1, datSResidual1t);
					djj+=1;
				}else
					fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
			}
			fclose(datSResidual1t);
			dll=sprintf(outdataname, "%s/BrawPMA_h%d.dat", OutputDir,i1); 
			fprintf(outf, "%s\n", outdataname);
			datSResidual1t=fopen(outdataname, "wb+");
			for(dkk=1, djj=1; dkk<=voxels; dkk++){
				if(flag[dkk]>=1)
				{
					if(waldpvalue[i1][djj]<=0.0)
					  mexPrintf("[%f] %f", waldpvalue[i1][djj], waldvalue[i1][djj]);
					fwrite(&waldpvalue[i1][djj], sizeof(float), 1, datSResidual1t);
			
					djj+=1;
				}else
					fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
			}
			fclose(datSResidual1t);
       		dll=sprintf(outdataname, "%s/BcorPMA_h%d.dat",OutputDir, i1); 
			fprintf(outf, "%s\n", outdataname);
			datSResidual1t=fopen(outdataname, "wb");
			for(dkk=1, djj=1; dkk<=voxels; dkk++){
				if(flag[dkk]>=1)
				{
					fwrite(&correctedwaldpvalue[i1][djj], sizeof(float), 1, datSResidual1t);
					djj+=1;
				}else
					fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
			}
			fclose(datSResidual1t);		
	   		
		}/* dmm */ 
	}/* end of if */

	else {
		for(dmm=1; dmm<=3; dmm++){     
			if(dmm==1)
				i1=1;
			if(dmm==2)
				i1=5;
			if(dmm==3)
				i1=hnum; 
			dll=sprintf(outdataname, "%s/BwaldMA_h%d.dat",OutputDir, i1);
			fprintf(outf, "%s\n", outdataname);
			datSResidual1t=fopen(outdataname, "wb+");
			for(dkk=1, djj=1; dkk<=voxels; dkk++){
				if(flag[dkk]>=1)
				{
					fwrite(&waldvalue[i1][djj], sizeof(float), 1, datSResidual1t);
					djj+=1;
				}else
					fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
			}
			fclose(datSResidual1t);
			dll=sprintf(outdataname, "%s/BrawPMA_h%d.dat", OutputDir,i1); 
			fprintf(outf, "%s\n", outdataname);
			datSResidual1t=fopen(outdataname, "wb+");
			for(dkk=1, djj=1; dkk<=voxels; dkk++){
				if(flag[dkk]>=1)
				{
					if(waldpvalue[i1][djj]<=0.0)
					  mexPrintf("[%f] %f", waldpvalue[i1][djj], waldvalue[i1][djj]);
					fwrite(&waldpvalue[i1][djj], sizeof(float), 1, datSResidual1t);
			
					djj+=1;
				}else
					fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
			}
			fclose(datSResidual1t);
       		dll=sprintf(outdataname, "%s/BcorPMA_h%d.dat",OutputDir, i1); 
			fprintf(outf, "%s\n", outdataname);
			datSResidual1t=fopen(outdataname, "wb+");
			for(dkk=1, djj=1; dkk<=voxels; dkk++){
				if(flag[dkk]>=1)
				{
					fwrite(&correctedwaldpvalue[i1][djj], sizeof(float), 1, datSResidual1t);
					djj+=1;
				}else
					fwrite(&temphongtu, sizeof(float), 1, datSResidual1t);
			}
			fclose(datSResidual1t);		
	   		
		}/* dmm */ 

	}/*end of else*/
	fclose(outf);
	delete[] outdataname;
	mxFree(OutputDir);


	free_vector(rr00, 0, DimYY*num_cov);
    free_ivector(flag,1,voxels);
	free_matrix(x0x, 0, TotalImg-1, 0, num_cov-1);
	free_dmatrix(designXX, 1, TotalImg*DimYY, 1, Ncov);
	free_dmatrix(allsurface, 1,NopointSurNew,1,TotalImg*DimYY);
	free_dvector(beta, 1,Ncov+2*DimYY);
	free_matrix(BrainTheta, 1, NopointSurNew, 1, Ncov+2*DimYY);	
	free_fmatrix3(BrainThetaVar, 1,NopointSurNew,1,Ncov,1,Ncov);
	free_dmatrix(linearbeta, 1,NopointSurNew,1,Ncov);
    free_ivector(Nrow, 1, Nsample);
    free_dvector(residualli, 1, totalobs);
	free_dmatrix(varGEE, 1, Ncov, 1,Ncov);
	free_fmatrix4(varMatrix, 1, NopointSurNew, 1, Nsample, 1, obsforonemax, 1, obsforonemax);
	free_dvector(waldvaluetemp, 1,NopointSurNew);
	free_fmatrix3(BrainThetaVarOld, 1,NopointSur,1,Ncov,1,Ncov);
    free_dmatrix(chivalue, 1,1,1,1);
	free_matrix(correctedwaldpvalue, 1,hnum,1,NopointSur);
	free_ivector(countbig, 1,hnum);
	free_ivector(countnumber, 1,NopointSur);
    free_dmatrix(diff, 1,Ncov,1,1);
	free_dmatrix(diff3, 1, Ncov, 1, 1);
	free_ivector(done, 1, NopointSur);
	free_dmatrix(h3beta, 1, NopointSur, 1, Ncov);
	free_matrix(holdbeta, 1, NopointSur, 1, Ncov);
	free_fmatrix3(INVvarGEE3, 1,NopointSur,1,Ncov,1,Ncov);
	free_imatrix(labelid, 1, NopointSur, 1, dllmax+1);
	free_matrix(sortwaldpvalue, 1,hnum, 1, NopointSur); 
	free_matrix(varGEE3D, 1, Ncov, 1, Ncov);
    free_matrix(waldpvalue, 1,hnum,1,NopointSur);
	free_dvector(waldpvalueold, 1,NopointSur);
	free_matrix(waldvalue, 1,hnum,1,NopointSur);
	free_dvector(waldvalueold, 1,NopointSur);
	free_dmatrix(weight, 1,NopointSur,1,dllmax+1);
	free_ivector(noneighs, 1, NopointSur);
	free_imatrix(neighID,1,NopointSur,1,300);
    free_matrix(distance,1,NopointSur,1,300); 

	mexPrintf("\n END OF MEXFUNCTION MAGEE.\n");
	dll=sprintf(logMesg, "END OF MEXFUNCTION MAGEE.");
	fprintf(runlog,"%s\n",logMesg);

    delete[] logMesg;
	fclose(runlog);

 	return;
			
}/* end */ 
