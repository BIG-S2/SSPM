
#ifndef PREDEF_HEADER
#define PREDEF_HEADER

#define NR_END 1
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

extern int  TotalImg, Ncov, NopointSur, Nsample,NcovVar, DimYY,DimXX, NcovDim, DimSPD,NOmaxTIMEs; 
extern unsigned int congrval,tausval;
extern int  htNoM,  htNoR,  *allDimTime;


void setSEED(int Y [],int SEED);
void init1(int seeds[]);
void nrerror(const char* error_text);
int *ivector(int nl, int nh);
unsigned long *lvector(int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float *vector(int nl, int nh);
double *dvector(int nl, int nh);
void chisqQuantile(float** chi2ref);  
float pythag(float a, float b);
void svdcmp(double **a,int m,int n,double *w, double **v);
void matrixmultiply (double **A, int row1,int column1,double **B, int row2,int column2,double **AB);
void matrixtranspose(double **A,int row,int column, double **B);
float ***fmatrix3(int nrl, int nrh, int ncl, int nch, int n3l, int n3h);
float ****fmatrix4(int nrl, int nrh, int ncl, int nch, int n3l, int n3h, int n4l, int n4h);
void readAnalyzeBody_char(char* fname, int dimz, int dimx, int dimy, float*** tmp);
void readAnalyzeDTI(char* fname, int dimz, int dimx, int dimy, int dimK, float*** tmp);
double DnewIVrank(double **DXXD, int *rankXX, int Nrow);
void matrixprint(double **AA,int Nrow,int Ncov);
void ABSmatrixNhalf(float **Corr, float **nHalfCorr, int NP, double powst);
void linearHongtu(double *beta,  double *residual, double **designXX,double **TXX, double *response);  
void CalMeanParaHongtu(double **designXX, int p, float ***varMatrix,   double *residual, double *beta, int *mi, int li,double *outputYYO, double *oldbeta);
void CalVarGEEHongtu(double **designXX, double *residual, int li, int *mi, int p, double *beta,   float ***varMatrix,  double **varGEEbeta);
void FinalCorrHongtu(float **Corr1, float **Corr2, float ***Final, int *mi, int li, int Nmax, float **ExactTime);
double AR1_dim(double x, double** xxYmatrix, int  N2total); 
double goldenHT(double ax, double bx, double cx, double (*f)(double, double**, int),  double tol, double *xmin, double **xxY, int N2total);
double EstimateAR_rho(double** timeDresidD, int* indxMi, int djj, int N2total);
void AR1TimeHongtu(double *residual, int *mi, int li, float **Corr, float **Corr2, float** ExactTime, int* indxMi, int* indxMi2);  
void CalVarWeightedGEECorrectedForVar(int current,double **outputYYO,double **designXX,int li,int *mi,float **beta,int q, double **varGEE,double *weight, int *labelid, int noneighs, float ****VarMatrix);
void GEEestimatesHongtu(double *beta,   double *residual, double **designXX, double **TXX, double *outputYYO, int li, int *mi, int p, int q, float ***varMatrix, float **ExactTime, int* indxMI, int* indxMI1);
float waldtest(float *BrainTheta, double **varGEE,float **RR,float *rr00,int noRow);
double gammln(double xx);
double  betacf(double a,double b,double x);
double betai(double a, double b, double x);
void CalMeanParaWeightGEECorrectedForVar (double **designXX, float *beta,int *mi, int li,double **outputYYO,double *weight,int *labelid,int noneighs,float ****VarMatrix);

void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_lvector(unsigned long *v, int nl, int nh);
void free_fmatrix3(float ***m, int nrl, int nrh, int ncl, int nch, int n3l, int n3h);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_fmatrix4(float ****m,int nrl, int nrh, int ncl, int nch, int n3l, int n3h, int n4l, int n4h);


int dmin( int a, int b);
int dmax( int a, int b);
void choldcYS (float **a, int n, float p[]);
void cholslYS(float **a,int n,float p[],float b[],float x[]);
void ivYS(float **W, int Nlow, int Nup);


#endif