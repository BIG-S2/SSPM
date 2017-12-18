# include <cstdio> 
# include <cstdlib>  
# include <cmath>    
# include <fstream>
# include <string.h>
#include "predef.h"
#include "mex.h"
#include "newdelete.h"


void setSEED(int Y [],int SEED){
	int X[]={21, 14, 49,0,1,2, 32, 22, 36, 23, 28,  3};
	int size=12;
	int i;
	unsigned int tmp=((SEED-1)%1024)*(int)pow(2.0,22.0);

	for(i=0;i<size;i++){
		Y[i]=X[i];
	}
	
	for(i=5;i>=3;i--){
		Y[i]=tmp/(int)pow(2.0,6.0*i);
		tmp-=Y[i]*(int)pow(2.0,6.0*i);
	}
}


//unsigned int congrval,tausval;
void init1(int seeds[]){
	int i;
	int size=6;
	unsigned int ctmp,ttmp;
	ctmp=0;
	ttmp=0;
	for(i=0;i<6;i++){
		ctmp+=seeds[i]*((unsigned int)pow(2.0,6.0*i));
		ttmp+=seeds[i+6]*((unsigned int)pow(2.0,6.0*i));
	}
	congrval=ctmp;
	tausval=ttmp;
}


void nrerror(const char* error_text)
{
//	void exit(int);

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	mexPrintf(error_text);
	//mexCallMATLAB(0,NULL,0,NULL,"pause");
	exit(1);
}


int *ivector(int nl, int nh)
{
	int *v;
   // v=(int *)mxCalloc(nh-nl+1,sizeof(int));
	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;

}


unsigned long *lvector(int nl, int nh)
/* allocate an unsigned int vector with subscript range v[nl..nh] */
{
	unsigned long *v;
    //v=(unsigned long *)mxCalloc(nh-nl+1+NR_END,sizeof(long));
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}


float **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	float **m;
     //m=(float **) mxCalloc(nrh-nrl+1,sizeof(float*));
	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		//m[i]=(float *) mxCalloc(nch-ncl+1,sizeof(float));
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}



int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;
    //m=(int **)mxCalloc(nrh-nrl+1,sizeof(int*));
	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		//m[i]=(int *)mxCalloc(nch-ncl+1,sizeof(int));
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



float *vector(int nl, int nh)
{
	float *v;
    v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}


double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}





void chisqQuantile(float** chi2ref)
{
	int dii, djj;
	float chi3ref[100][9] =
	{ 0.0158, 0.0642, 0.1485, 0.2750, 0.4549, 0.7083, 1.0742, 1.6424, 2.7055,
	0.2107, 0.4463, 0.7133, 1.0217, 1.3863, 1.8326, 2.4079, 3.2189, 4.6052,
	0.5844, 1.0052, 1.4237, 1.8692, 2.3660, 2.9462, 3.6649, 4.6416, 6.2514,
	1.0636, 1.6488, 2.1947, 2.7528, 3.3567, 4.0446, 4.8784, 5.9886, 7.7794,
	1.6103, 2.3425, 2.9999, 3.6555, 4.3515, 5.1319, 6.0644, 7.2893, 9.2364,
	2.2041, 3.0701, 3.8276, 4.5702, 5.3481, 6.2108, 7.2311, 8.5581, 10.6446,
	2.8331, 3.8223, 4.6713, 5.4932, 6.3458, 7.2832, 8.3834, 9.8032, 12.0170,
	3.4895, 4.5936, 5.5274, 6.4226, 7.3441, 8.3505, 9.5245, 11.0301, 13.3616,
	4.1682, 5.3801, 6.3933, 7.3570, 8.3428, 9.4136, 10.6564, 12.2421, 14.6837,
	4.8652, 6.1791, 7.2672, 8.2955, 9.3418, 10.4732, 11.7807, 13.4420, 15.9872,
	5.5778, 6.9887, 8.1479, 9.2373, 10.3410, 11.5298, 12.8987, 14.6314, 17.2750,
	6.3038, 7.8073, 9.0343, 10.1820, 11.3403, 12.5838, 14.0111, 15.8120, 18.5493,
	7.0415, 8.6339, 9.9257, 11.1291, 12.3398, 13.6356, 15.1187, 16.9848, 19.8119,
	7.7895, 9.4673, 10.8215, 12.0785, 13.3393, 14.6853, 16.2221, 18.1508, 21.0641,
	8.5468, 10.3070, 11.7212, 13.0297, 14.3389, 15.7332, 17.3217, 19.3107, 22.3071,
	9.3122, 11.1521, 12.6243, 13.9827, 15.3385, 16.7795, 18.4179, 20.4651, 23.5418,
	10.0852, 12.0023, 13.5307, 14.9373, 16.3382, 17.8244, 19.5110, 21.6146, 24.7690,
	10.8649, 12.8570, 14.4399, 15.8932, 17.3379, 18.8679, 20.6014, 22.7595, 25.9894,
	11.6509, 13.7158, 15.3517, 16.8504, 18.3377, 19.9102, 21.6891, 23.9004, 27.2036,
	12.4426, 14.5784, 16.2659, 17.8088, 19.3374, 20.9514, 22.7745, 25.0375, 28.4120,
	13.2396, 15.4446, 17.1823, 18.7683, 20.3372, 21.9915, 23.8578, 26.1711, 29.6151,
	14.0415, 16.3140, 18.1007, 19.7288, 21.3370, 23.0307, 24.9390, 27.3015, 30.8133,
	14.8480, 17.1865, 19.0211, 20.6902, 22.3369, 24.0689, 26.0184, 28.4288, 32.0069,
	15.6587, 18.0618, 19.9432, 21.6525, 23.3367, 25.1063, 27.0960, 29.5533, 33.1962,
	16.4734, 18.9398, 20.8670, 22.6156, 24.3366, 26.1430, 28.1719, 30.6752, 34.3816,
	17.2919, 19.8202, 21.7924, 23.5794, 25.3365, 27.1789, 29.2463, 31.7946, 35.5632,
	18.1139, 20.7030, 22.7192, 24.5440, 26.3363, 28.2141, 30.3193, 32.9117, 36.7412,
	18.9392, 21.5880, 23.6475, 25.5093, 27.3362, 29.2486, 31.3909, 34.0266, 37.9159,
	19.7677, 22.4751, 24.5770, 26.4751, 28.3361, 30.2825, 32.4612, 35.1394, 39.0875,
	20.5992, 23.3641, 25.5078, 27.4416, 29.3360, 31.3159, 33.5302, 36.2502, 40.2560,
	21.4336, 24.2551, 26.4397, 28.4087, 30.3359, 32.3486, 34.5981, 37.3591, 41.4217,
	22.2706, 25.1478, 27.3728, 29.3763, 31.3359, 33.3809, 35.6649, 38.4663, 42.5847,
	23.1102, 26.0422, 28.3069, 30.3444, 32.3358, 34.4126, 36.7307, 39.5718, 43.7452,
	23.9523, 26.9383, 29.2421, 31.3130, 33.3357, 35.4438, 37.7954, 40.6756, 44.9032,
	24.7967, 27.8359, 30.1782, 32.2821, 34.3356, 36.4746, 38.8591, 41.7780, 46.0588,
	25.6433, 28.7350, 31.1152, 33.2517, 35.3356, 37.5049, 39.9220, 42.8788, 47.2122,
	26.4921, 29.6355, 32.0532, 34.2216, 36.3355, 38.5348, 40.9839, 43.9782, 48.3634,
	27.3430, 30.5373, 32.9919, 35.1920, 37.3355, 39.5643, 42.0451, 45.0763, 49.5126,
	28.1958, 31.4405, 33.9315, 36.1628, 38.3354, 40.5935, 43.1053, 46.1730, 50.6598,
	29.0505, 32.3450, 34.8719, 37.1340, 39.3353, 41.6222, 44.1649, 47.2685, 51.8051,
	29.9071, 33.2506, 35.8131, 38.1055, 40.3353, 42.6506, 45.2236, 48.3628, 52.9485,
	30.7654, 34.1574, 36.7550, 39.0774, 41.3352, 43.6786, 46.2817, 49.4560, 54.0902,
	31.6255, 35.0653, 37.6975, 40.0496, 42.3352, 44.7063, 47.3390, 50.5480, 55.2302,
	32.4871, 35.9743, 38.6408, 41.0222, 43.3352, 45.7336, 48.3957, 51.6389, 56.3685,
	33.3504, 36.8844, 39.5847, 41.9950, 44.3351, 46.7607, 49.4517, 52.7288, 57.5053,
	34.2152, 37.7955, 40.5292, 42.9682, 45.3351, 47.7874, 50.5071, 53.8177, 58.6405,
	35.0814, 38.7075, 41.4744, 43.9417, 46.3350, 48.8139, 51.5619, 54.9056, 59.7743,
	35.9491, 39.6205, 42.4201, 44.9154, 47.3350, 49.8401, 52.6161, 55.9926, 60.9066,
	36.8182, 40.5344, 43.3664, 45.8895, 48.3350, 50.8660, 53.6697, 57.0786, 62.0375,
	37.6886, 41.4492, 44.3133, 46.8638, 49.3349, 51.8916, 54.7228, 58.1638, 63.1671,
	38.5604, 42.3649, 45.2607, 47.8383, 50.3349, 52.9170, 55.7753, 59.2481, 64.2954,
	39.4334, 43.2814, 46.2086, 48.8132, 51.3349, 53.9421, 56.8274, 60.3316, 65.4224,
	40.3076, 44.1987, 47.1571, 49.7882, 52.3348, 54.9670, 57.8789, 61.4142, 66.5482,
	41.1830, 45.1167, 48.1060, 50.7635, 53.3348, 55.9916, 58.9299, 62.4961, 67.6728,
	42.0596, 46.0356, 49.0554, 51.7391, 54.3348, 57.0160, 59.9805, 63.5772, 68.7962,
	42.9373, 46.9552, 50.0053, 52.7148, 55.3348, 58.0402, 61.0305, 64.6576, 69.9185,
	43.8161, 47.8755, 50.9556, 53.6908, 56.3347, 59.0642, 62.0802, 65.7373, 71.0397,
	44.6960, 48.7965, 51.9063, 54.6670, 57.3347, 60.0879, 63.1294, 66.8162, 72.1598,
	45.5770, 49.7182, 52.8575, 55.6434, 58.3347, 61.1115, 64.1782, 67.8945, 73.2789,
	46.4589, 50.6406, 53.8091, 56.6200, 59.3347, 62.1348, 65.2265, 68.9721, 74.3970,
	47.3418, 51.5636, 54.7611, 57.5968, 60.3346, 63.1580, 66.2745, 70.0490, 75.5141,
	48.2257, 52.4873, 55.7135, 58.5738, 61.3346, 64.1810, 67.3220, 71.1253, 76.6302,
	49.1105, 53.4116, 56.6663, 59.5510, 62.3346, 65.2037, 68.3692, 72.2010, 77.7454,
	49.9963, 54.3365, 57.6195, 60.5283, 63.3346, 66.2263, 69.4160, 73.2761, 78.8596,
	50.8829, 55.2620, 58.5731, 61.5059, 64.3346, 67.2488, 70.4624, 74.3506, 79.9730,
	51.7705, 56.1880, 59.5270, 62.4836, 65.3345, 68.2710, 71.5085, 75.4245, 81.0855,
	52.6588, 57.1147, 60.4812, 63.4615, 66.3345, 69.2931, 72.5542, 76.4978, 82.1971,
	53.5481, 58.0418, 61.4358, 64.4395, 67.3345, 70.3150, 73.5995, 77.5707, 83.3079,
	54.4381, 58.9696, 62.3908, 65.4177, 68.3345, 71.3368, 74.6446, 78.6429, 84.4179,
	55.3289, 59.8978, 63.3460, 66.3961, 69.3345, 72.3583, 75.6893, 79.7146, 85.5270,
	56.2206, 60.8266, 64.3016, 67.3746, 70.3345, 73.3798, 76.7337, 80.7859, 86.6354,
	57.1129, 61.7558, 65.2575, 68.3533, 71.3344, 74.4011, 77.7777, 81.8566, 87.7430,
	58.0061, 62.6856, 66.2137, 69.3322, 72.3344, 75.4222, 78.8215, 82.9268, 88.8499,
	58.9000, 63.6158, 67.1702, 70.3111, 73.3344, 76.4432, 79.8650, 83.9965, 89.9560,
	59.7946, 64.5466, 68.1271, 71.2903, 74.3344, 77.4640, 80.9081, 85.0658, 91.0615,
	60.6899, 65.4777, 69.0842, 72.2695, 75.3344, 78.4848, 81.9510, 86.1346, 92.1662,
	61.5858, 66.4094, 70.0415, 73.2489, 76.3344, 79.5053, 82.9936, 87.2030, 93.2702,
	62.4825, 67.3415, 70.9992, 74.2285, 77.3344, 80.5258, 84.0359, 88.2709, 94.3735,
	63.3799, 68.2740, 71.9571, 75.2081, 78.3343, 81.5461, 85.0779, 89.3383, 95.4762,
	64.2778, 69.2069, 72.9153, 76.1879, 79.3343, 82.5663, 86.1197, 90.4053, 96.5782,
	65.1765, 70.1403, 73.8738, 77.1679, 80.3343, 83.5863, 87.1612, 91.4720, 97.6796,
	66.0757, 71.0741, 74.8325, 78.1479, 81.3343, 84.6062, 88.2025, 92.5382, 98.7803,
	66.9756, 72.0083, 75.7915, 79.1281, 82.3343, 85.6260, 89.2435, 93.6039, 99.8805,
	67.8761, 72.9429, 76.7507, 80.1084, 83.3343, 86.6457, 90.2842, 94.6693, 100.9800,
	68.7772, 73.8779, 77.7102, 81.0888, 84.3343, 87.6653, 91.3247, 95.7343, 102.0789,
	69.6788, 74.8132, 78.6699, 82.0693, 85.3343, 88.6847, 92.3650, 96.7990, 103.1773,
	70.5810, 75.7490, 79.6299, 83.0500, 86.3342, 89.7041, 93.4050, 97.8632, 104.2750,
	71.4838, 76.6851, 80.5901, 84.0307, 87.3342, 90.7233, 94.4448, 98.9271, 105.3722,
	72.3872, 77.6216, 81.5505, 85.0116, 88.3342, 91.7424, 95.4844, 99.9906, 106.4689,
	73.2911, 78.5584, 82.5111, 85.9925, 89.3342, 92.7614, 96.5238, 101.0537, 107.5650,
	74.1955, 79.4956, 83.4719, 86.9736, 90.3342, 93.7803, 97.5629, 102.1165, 108.6606,
	75.1005, 80.4332, 84.4330, 87.9548, 91.3342, 94.7991, 98.6018, 103.1790, 109.7556,
	76.0060, 81.3711, 85.3943, 88.9361, 92.3342, 95.8178, 99.6405, 104.2411, 110.8502,
	76.9119, 82.3093, 86.3558, 89.9175, 93.3342, 96.8364, 100.6790, 105.3028, 111.9442,
	77.8184, 83.2478, 87.3175, 90.8990, 94.3342, 97.8549, 101.7173, 106.3643, 113.0377,
	78.7254, 84.1867, 88.2794, 91.8806, 95.3342, 98.8733, 102.7554, 107.4254, 114.1307,
	79.6329, 85.1259, 89.2415, 92.8622, 96.3342, 99.8916, 103.7933, 108.4862, 115.2232,
	80.5408, 86.0654, 90.2038, 93.8440, 97.3341, 100.9098, 104.8310, 109.5467, 116.3153,
	81.4493, 87.0052, 91.1663, 94.8259, 98.3341, 101.9279, 105.8685, 110.6068, 117.4069,
	82.3581, 87.9453, 92.1289, 95.8078, 99.3341, 102.9459, 106.9058, 111.6667, 118.4980
	};

	for (dii = 1; dii <= 100; dii++)
	for (djj = 1; djj <= 9; djj++)
		chi2ref[dii][djj] = chi3ref[dii - 1][djj - 1];

}/* end  */

	
	
float pythag(float a, float b)
{
  float absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) 
	  return absa*sqrt(1.0+SQR(absb/absa));
  else 
	  return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}/* end */



void svdcmp(double **a,int m,int n,double *w, double **v)
{
	int flag,i,its,j,jj,k,l,nm;
	double c,f,h,s,x,y,z;
	double anorm=0.0,g=0.0,scale=0.0;
	double *rv1;
 
	if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
	rv1=dvector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=(float)(1.0/g);
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=60;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=pythag(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
//			if (its == 30) nrerror("No convergence in 30 SVDCMP iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}



void matrixmultiply (double **A, int row1,int column1,double **B, int row2,int column2,double **AB)
{
	int i,j,m;
//	float **AB;
//	AB=dmatrix(1,row1,1,column2);
	for (i=1;i<=row1;i++)
		for (j=1;j<=column2;j++)
			AB[i][j]=0;

	for (i=1;i<=row1;i++)
  for (j=1;j<=column2;j++)
	  for (m=1;m<=column1;m++){
      AB[i][j]+=A[i][m]*B[m][j];
//printf("\n A=%14f B=%14f AB=%14f",A[i][m],B[m][j],AB[i][j]);
	  }
   //return(AB);
   //free_dmatrix(AB,1,row1,1,column2);
}
void matrixtranspose(double **A,int row,int column, double **B)
{ 
	int i,j;
//	float **B;
//	B=dmatrix(1,column,1,row);
     for (i=1;i<=column;i++)
		for (j=1;j<=row;j++)
			B[i][j]=0;

   for  (i=1;i<=row;i++)
	   for (j=1;j<=column;j++){
           B[j][i]=A[i][j];
//	printf("B=%14lf",B[i][j]);
	   }
//		return(B);
//	free_dmatrix(B,1,column,1,row);
}

// mask_tmp = fmatrix3(0, dimZ-1, 0, dimX-1, 0, dimY-1);
float ***fmatrix3(int nrl, int nrh, int ncl, int nch, int n3l, int n3h)
{
  int i, j;
  float ***m;
  m=(float ***) malloc((unsigned) (nrh-nrl+1)*sizeof(float**));
  if (!m) nrerror("error in fmatrix3");
  m -= nrl;
  for(i=nrl; i<=nrh;i++) {
    m[i]=(float **) malloc((unsigned) (nch-ncl+1)*sizeof(float*));
    if (!m[i]) nrerror("error in dmatrix3");
    m[i] -= ncl;
    for(j=ncl; j<=nch; j++){
      m[i][j]=(float *) malloc((unsigned) (n3h-n3l+1)*sizeof(float));
      if (!m[i][j]) nrerror("error in dmatrix3");
	  m[i][j] -= n3l;
    }      
  }
  return m;
}


float ****fmatrix4(int nrl, int nrh, int ncl, int nch, int n3l, int n3h, int n4l, int n4h)
{
  int i, j, k;
  float ****m;

  m=(float ****) malloc((unsigned) (nrh-nrl+1)*sizeof(float***));
  if (!m) nrerror("error in fmatrix4");
  m -= nrl;
  for(i=nrl; i<=nrh;i++) {
    m[i]=(float ***) malloc((unsigned) (nch-ncl+1)*sizeof(float**));
    if (!m[i]) nrerror("error in fmatrix4");
    m[i] -= ncl;
    for(j=ncl; j<=nch; j++){
      m[i][j]=(float **) malloc((unsigned) (n3h-n3l+1)*sizeof(float*));
      if (!m[i][j]) nrerror("error in fmatrix4");
	  m[i][j] -= n3l;
	  for(k=n3l; k<=n3h; k++){
         m[i][j][k]=(float *) malloc((unsigned) (n4h-n4l+1)*sizeof(float));
         if (!m[i][j][k]) nrerror("error in fmatrix4");
	     m[i][j][k] -= n4l;
	  }
    }      
  }
  return m;
}



void free_vector(float *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_ivector(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}




void free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}


void free_lvector(unsigned long *v, int nl, int nh)
/* free an unsigned int vector allocated with lvector() */
{
	free((char*)(v+nl-NR_END));
	//free((FREE_ARG)(v+nl-NR_END));
}



void free_fmatrix3(float ***m, int nrl, int nrh, int ncl, int nch, int n3l, int n3h) 
{
	int i,j;
	for(i=nrh;i>=nrl;i--)
   	   for(j=nch;j>=ncl;j--)
		free((char*) (m[i][j]+n3l));
		for(i=nrh;i>=nrl;i--) 
		free((char*) (m[i]+ncl));
		free((char*) (m+nrl));
}

void free_fmatrix4(float ****m,int nrl, int nrh, int ncl, int nch, int n3l, int n3h, int n4l, int n4h) 
{
	int i,j, k;
	for(i=nrh;i>=nrl;i--)
   	   for(j=nch;j>=ncl;j--)
		   for(k=n3h; k>=n3l; k--)
		free((char*) (m[i][j][k]+n4l));
	for(i=nrh;i>=nrl;i--)
   	   for(j=nch;j>=ncl;j--)
		free((char*) (m[i][j]+n3l));
	for(i=nrh;i>=nrl;i--) 
		free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) //mxFree(m[i]+ncl);
	free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) 
	free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


//readAnalyzeBody_char((char*)maskTempo, dimZ, dimX, dimY, mask_tmp);
void readAnalyzeBody_char(char* fname, int dimz, int dimx, int dimy, float*** tmp)
{
int i, j, k, totalBYTE;
FILE* fidread = fopen(fname, "rb");
if(fidread==NULL){
   mexPrintf("Cannot open the save file\n");
   mexPrintf("filename is:%s\n",fname);
   mexCallMATLAB(0,NULL,0,NULL,"pause");
   exit(1);//////////
    }
//totalBYTE=dimx*dimy*dimz*4;///////////////////////////////////////////////////
totalBYTE = 256*256*256*24;
char* buffer = new char[totalBYTE];
int read = fread(buffer, sizeof(char), totalBYTE, fidread);
fclose(fidread);
delete[] buffer;
int byteSrc = read/dimz/dimx/dimy; // # of bytes for each data value

fidread = fopen(fname, "rb");
switch(byteSrc){

case 1:{ 
unsigned char* img= new unsigned char[dimx];
for(k=0; k<dimz; k++)
for(j=0; j<dimy; j++)
{
fread(img, sizeof(unsigned char), dimx, fidread);
for(i=0; i<dimx; i++)
tmp[k][i][j] = img[i]; 
}/* for */ 
delete[] img;
      }
break;

case 2:{   

unsigned short* img= new unsigned short[dimx];
for(k=0; k<dimz; k++)
for (j = 0; j<dimy; j++)
{
	fread(img, sizeof(unsigned short), dimx, fidread);
	for (i = 0; i<dimx; i++)
		tmp[k][i][j] = img[i];
}/* for */ 
delete[] img;
	   }
	   break;

case 4:{  

float* img =new float[dimx];
for(k=0; k<dimz; k++)
for (j = 0; j<dimy; j++)
{
	fread(img, sizeof(float), dimx, fidread);
	for (i = 0; i<dimx; i++)
		tmp[k][i][j] = img[i];
}/* for */ 
delete[] img;
	   } 
	   break;

default:{
mexPrintf("unsupported data size with  %d bytes long.\n", byteSrc);
mexCallMATLAB(0,NULL,0,NULL,"pause");
exit(1);
}
break;
}/* switch */  
fclose(fidread);
}


//readAnalyzeDTI(name_input_img, dimZ, dimX, dimY, 1, tensors);
void readAnalyzeDTI(char* fname, int dimz, int dimx, int dimy, int dimK, float*** tmp)
{
int i, j, k, kk, totalBYTE, indx;
FILE* fidread = fopen(fname, "rb");
if(fidread==NULL){
   mexPrintf("Cannot open the save file\n");
   mexPrintf("filename is:%s\n",fname);
   mexCallMATLAB(0,NULL,0,NULL,"pause");
   exit(1);//////////
    }

//totalBYTE=dimx*dimy*dimz*dimK*4;////////////////////////////////
totalBYTE = 256*256*256*24;
char* buffer = new char[totalBYTE];
int read = fread(buffer, sizeof(char), totalBYTE, fidread);
fclose(fidread);
delete[] buffer;
int byteSrc = read/dimz/dimx/dimy/dimK; // # of bytes for each data value

fidread = fopen(fname, "rb");
int voxNUM=dimx*dimy*dimz*dimK; 

switch(byteSrc){

case 1:{ 
unsigned char* img = new unsigned char[voxNUM];
fread(img, sizeof(unsigned char), voxNUM, fidread);
indx=0; 
for(k=0; k<dimz; k++)
for (j = 0; j<dimy; j++)
for (i = 0; i<dimx; i++)
{
for(kk=0; kk<dimK; kk++, indx++) {
tmp[k*dimK+kk][i][j] = img[indx];
//   printf("\n tmp=%f",tmp[k*dimK+kk][i][j]);
}
}
delete[] img;
	   } 
	   break;


case 2:{ 
unsigned short* img = new unsigned short[voxNUM];
fread(img, sizeof(unsigned short), voxNUM, fidread);
indx=0; 
for(k=0; k<dimz; k++)
for (j = 0; j<dimy; j++)
for (i = 0; i<dimx; i++)
{
for(kk=0; kk<dimK; kk++, indx++) {
tmp[k*dimK+kk][i][j] =img[indx];
//   printf("\n tmp=%f",tmp[k*dimK+kk][i][j]);
}
}
delete[] img;
} 
	   break;


case 4:{  
float* img = new float[voxNUM];
fread(img, sizeof(float), voxNUM, fidread);
indx=0; 
for(k=0; k<dimz; k++)
for (j = 0; j<dimy; j++)
for (i = 0; i<dimx; i++)
{
for(kk=0; kk<dimK; kk++, indx++) {
tmp[k*dimK+kk][i][j] = img[indx];
//   printf("\n tmp=%f",tmp[k*dimK+kk][i][j]);
}
}
delete[] img;
} 
	   break;

default:{
mexPrintf("unsupported data size with  %d bytes long.\n", byteSrc);
mexCallMATLAB(0,NULL,0,NULL,"pause");
exit(1);
 }
break;

}/* switch */ 
fclose(fidread);	  
}



double DnewIVrank(double **DXXD, int *rankXX,   int Nrow)
{
 	int dkk, djj, dii; 
	double determ; 
	double *WW, **VV, **UU, **PP1, **PP0, tempt, sumEigen; 
	
	WW=dvector(1, Nrow+1);
	UU=dmatrix(1, Nrow+1, 1, Nrow+1); 
	VV=dmatrix(1, Nrow+1, 1, Nrow+1); 
	PP1=dmatrix(1, Nrow+1, 1, Nrow+1);
	PP0=dmatrix(1, Nrow+1, 1, Nrow+1); 
 
	for(dkk=1; dkk<=Nrow; dkk++)
		for(djj=1; djj<=Nrow; djj++){
			UU[dkk][djj]=DXXD[dkk][djj]*1.0;
			PP1[dkk][djj]=0.0;
		    PP0[dkk][djj]=0.0; 
		  }
    svdcmp(UU, Nrow, Nrow, WW, VV); 
	tempt=1.0;
	sumEigen=0.0; 
	for(dii=1; dii<=Nrow; dii++){
  	   tempt=tempt*WW[dii]; 
       sumEigen+=fabs(WW[dii]);
	}
	sumEigen=sumEigen/(Nrow*1.0); 
	determ=(float)tempt; 
 
    rankXX[1]=0; 
    for(dii=1; dii<=Nrow; dii++)
 	    if(fabs(WW[dii])>0.00000000000000000000000001*sumEigen&&fabs(WW[dii])>0){ 
	       WW[dii]=1.0/WW[dii]; 
           rankXX[1]++;   
		}else 
	       WW[dii]=0.0; 
		
    for(dii=1; dii<=Nrow; dii++) 
 	  for(djj=1; djj<=dii; djj++){
        DXXD[dii][djj]=0.0;
		tempt=0.0; 
	    for(dkk=1; dkk<=Nrow; dkk++)
           tempt+=VV[dii][dkk]*WW[dkk]*UU[djj][dkk]; 
		DXXD[dii][djj]=tempt; 
		DXXD[djj][dii]=tempt;  
	  }

	 free_dvector(WW, 1, Nrow+1);
	 free_dmatrix(UU, 1, Nrow+1, 1, Nrow+1); 
	 free_dmatrix(VV, 1, Nrow+1, 1, Nrow+1); 
	 free_dmatrix(PP1, 1, Nrow+1, 1, Nrow+1);
	 free_dmatrix(PP0, 1, Nrow+1, 1, Nrow+1); 	
     return(determ); 

}/* end */ 


void matrixprint(double **AA,int Nrow,int Ncov)
{

	int dii, djj;

	mexPrintf("\n");
	for (dii=1;dii<=Nrow;dii++){
		for(djj=1;djj<=Ncov;djj++)
		{
			mexPrintf("%f,",AA[dii][djj]);
		}
		mexPrintf("\n");
}

}

void ABSmatrixNhalf(float **Corr, float **nHalfCorr, int NP, double powst)
{ 
	double  **PP0,**PP1,**UU,**VV,*WW, tempt, LAMmin, LAMmax;
	int dii, dkk, djj, rankXX, NoP0;
	
 	PP0=dmatrix(1, NP, 1, NP);
 	PP1=dmatrix(1, NP, 1, NP);
 	UU=dmatrix(1, NP, 1, NP);
 	VV=dmatrix(1, NP, 1, NP);
 	WW=dvector(1, NP);
	
	for(dii=1; dii<=NP; dii++)
		for(djj=1; djj<=NP; djj++){
			UU[dii][djj]=Corr[dii][djj];
			PP1[dii][djj]=0.0;
			PP0[dii][djj]=0.0; 
		}
	
    svdcmp(UU, NP, NP, WW, VV); 
	rankXX=0; 
	NoP0=0; 
 
	LAMmin=0.0;
	LAMmax=0.0;
	for(dii=1; dii<=NP; dii++) 
	    if(fabs(WW[dii])>0.0){
		   if(fabs(WW[dii])>LAMmax) 
			   LAMmax=fabs(WW[dii]); 
		}
	LAMmin=LAMmax; 
	for(dii=1; dii<=NP; dii++) 
	    if(fabs(WW[dii])>0){
			if(fabs(WW[dii])<LAMmin) 
				LAMmin=fabs(WW[dii]); 
		}
	float precision=log10(LAMmax)-log10(LAMmin); 
	if(precision>6) 
		for(dii=1; dii<=NP; dii++) 
			if(fabs(WW[dii])>0){	
				WW[dii]=WW[dii]+LAMmax*0.000001;
			}	
	
    for(dii=1; dii<=NP; dii++){
		if(fabs(WW[dii])>0.000000001){ 
			WW[dii]=pow(fabs(WW[dii]), powst); 
			rankXX++;   
		}else{
			WW[dii]=0.0; 
			NoP0++; 
			for(djj=1; djj<=NP; djj++) 
				PP0[djj][NoP0]=VV[djj][dii];      
		}   
	}
	
	
    for(dii=1; dii<=NP; dii++) 
		for(djj=1; djj<=dii; djj++){
			tempt=0.0; 
			for(dkk=1; dkk<=NP; dkk++)
				tempt+=VV[dii][dkk]*WW[dkk]*VV[djj][dkk];
			nHalfCorr[dii][djj]=tempt; 
			nHalfCorr[djj][dii]=tempt; 
    	}
	free_dmatrix(PP0, 1, NP, 1, NP);
 	free_dmatrix(PP1, 1, NP, 1, NP);
 	free_dmatrix(UU, 1, NP, 1, NP);
 	free_dmatrix(VV, 1, NP, 1, NP);
 	free_dvector(WW, 1, NP);
		
}/* end */ 


void linearHongtu(double *beta,  double *residual, double **designXX,double **TXX, double *response)  
{
	int dii, djj, dkk;
	double  *XTY; 
    double SSE;

	XTY=dvector(1, Ncov+1); 

	
//	for(dkk=1; dkk<=TotalImg*DimSPD; dkk++)
//		printf("\n response[dkk]=%f",response[dkk]);

   
	for(dii=1; dii<=Ncov; dii++){
		XTY[dii]=0.0;
		for(dkk=1; dkk<=TotalImg*DimSPD; dkk++)
            XTY[dii]+=designXX[dkk][dii]*response[dkk];
	}


	for(dii=1; dii<=Ncov; dii++){
	    beta[dii]=0.0;
		for(djj=1; djj<=Ncov; djj++)
            beta[dii]+=TXX[dii][djj]*XTY[djj];
	//	printf("\nbeta[dii]=%f",beta[dii]);
	}

	SSE=0.0;
	for(dii=1; dii<=TotalImg*DimSPD; dii++)
		SSE+=response[dii]*response[dii];
    for(dii=1; dii<=Ncov; dii++)
	    SSE-=beta[dii]*XTY[dii]; 


	for(dii=1; dii<=TotalImg*DimSPD; dii++){
        residual[dii]=response[dii];
		for(dkk=1; dkk<=Ncov; dkk++)
			residual[dii]-=designXX[dii][dkk]*beta[dkk]; 
		
	}
    beta[Ncov+1]=SSE; 
    free_dvector(XTY, 1, Ncov+1);
  	

}/* end */ 



void AR1TimeHongtu(double *residual, int *mi, int li, float **Corr, float **Corr2, float** ExactTime, int* indxMi, int* indxMi2)  
{
	int dii, djj, dkk, dll, dmm, dnn, dtt,  position1, position2, N1total;
	int N2total=0; 
	double sum, alpha,phi;
	float **nHalfCorr, **Corr3; 
	double **CORRresidual; 
	double **timeDresidD; 

    for(dii=1; dii<=Nsample; dii++) 
	   N2total+=mi[dii]*(mi[dii]-1);
    N2total=(int)(N2total/2.0); 
   
	Corr3=matrix(1, li+1, 1, li+1); 
	CORRresidual=dmatrix(1, TotalImg, 1, DimSPD); 
	nHalfCorr=matrix(1, li+1, 1, li+1); 
    timeDresidD=dmatrix(1, N2total, 1, DimSPD+1);    

    for(dii=1; dii<=li; dii++)
		for(djj=1; djj<=li; djj++){
			Corr[dii][djj]=0.0; 
    	    Corr2[dii][djj]=0.0; 
		}
    
	for(dii=1; dii<=Nsample; dii++)
       for(djj=1; djj<=mi[dii]; djj++){
		   for(dkk=1; dkk<=li; dkk++)
			   for(dll=1; dll<=li; dll++){ 
				   position1=indxMi[dii-1]*li+li*(djj-1)+dkk;
                   position2=indxMi[dii-1]*li+li*(djj-1)+dll;
	 	           Corr[dkk][dll]+=residual[position1]*residual[position2];
			   } 
	 }
			   
	for(djj=1; djj<=li; djj++)
		for(dkk=1; dkk<=li; dkk++)
				Corr[djj][dkk]=Corr[djj][dkk]/(TotalImg*1.0-Ncov);
 
 
	ABSmatrixNhalf(Corr, nHalfCorr, li, -0.5); 


   /* new AR1 process Hongtu Dec 6 */ 
    dmm=1; 
    for(dii=1; dii<=Nsample; dii++){
	   for(djj=1; djj<=mi[dii]; djj++){ 
	   	   for(dkk=1; dkk<=li; dkk++){ 
			   CORRresidual[dmm][dkk]=0.0; 
			   for(dtt=1; dtt<=li; dtt++){
		          position1=indxMi[dii-1]*li+li*(djj-1)+dtt;
				  CORRresidual[dmm][dkk]+=nHalfCorr[dkk][dtt]*residual[position1]; 
			   }
		   }/* dkk */ 
            dmm++; 
	     }/* djj */  
	 }/* for */ 
   

   dmm=1; 
   for(dii=1; dii<=Nsample; dii++){
	   for(djj=1; djj<=mi[dii]; djj++)
	      for(dll=djj+1; dll<=mi[dii]; dll++){ 
		   position1=indxMi[dii-1]+djj;
           position2=indxMi[dii-1]+dll;
           timeDresidD[dmm][1]=fabs(ExactTime[dii][djj]-ExactTime[dii][dll]); 
		   for(dkk=1; dkk<=li; dkk++) 
              timeDresidD[dmm][dkk+1]=CORRresidual[position1][dkk]*CORRresidual[position2][dkk];
		   dmm++; 
	    }/* djj */ 
   }/* dii */ 

 

   if(N2total>0){
	   for(djj=1; djj<=li; djj++){
          alpha=EstimateAR_rho(timeDresidD, indxMi, djj, N2total); 
          Corr2[djj][djj]=(float)alpha; 
	   }
   }else{
     for(djj=1; djj<=li; djj++)
       Corr2[djj][djj]=0.0; 
   }     

	 free_matrix(Corr3, 1, li+1, 1, li+1); 
	 free_dmatrix(CORRresidual, 1, TotalImg, 1, DimSPD); 
	 free_matrix(nHalfCorr, 1, li+1, 1, li+1); 
     free_dmatrix(timeDresidD, 1, N2total, 1, DimSPD+1);    

}/* end */ 


  void	CalVarGEEHongtu(double **designXX, double *residual, int li, int *mi, int p, double *beta,   float ***varMatrix,  double **varGEEbeta)
  {

 	double  **inter1, **inter2, **inter3, **inter3T, **inter4, **inter5, **inter6, **inter7, **sum1, **sum2, **sum1T;
	double  **term1, **term2, **term3, **term4, **corrMatrixI, **designXXI;
	double  **designXXIT, **varMatrixI, **residualD, **residualDI, **residualDIT, **PP0, **PP1, **UU, *WW, **VV, *ZZxi, **covYI, **residualxiI, **residualxiIT;
	double  **ParVarI, **ParVarIT, **sum3, **sum4,**sum3T, **inter8;
    int     dii,djj,dkk,dll,dmm, Ntotal, position, Npar, mimax, Nrowmax, *Nrow, Nrowtotal, rankXX, NoP0; 
	float   ***varMatrix2, **sumtemp1, **sumtemp2;
	mimax=0;
	Nrowmax=0;
	
	Nrow=ivector(1,Nsample);
	Nrowtotal=0;


	for (dii=1;dii<=Nsample;dii++){
		Nrow[dii]=mi[dii]*li;
		if (mi[dii]>mimax) mimax=mi[dii];
	    Nrowtotal+=Nrow[dii];
	}

     Nrowmax=mimax*li;
 
	covYI=dmatrix(1,Nrowmax,1,Nrowmax);
	designXXI=dmatrix(1,Nrowmax,1,Ncov);
	designXXIT=dmatrix(1,Ncov,1,Nrowmax);
 	inter3=dmatrix(1,Ncov,1,Nrowmax);
	inter3T=dmatrix(1,Nrowmax,1,Ncov);
	inter4=dmatrix(1,Ncov,1,Nrowmax);
	inter5=dmatrix(1,Ncov,1,Ncov);
  	residualD=dmatrix(1,Nrowtotal,1,1);
	residualDI=dmatrix(1,Nrowmax,1,1);
	residualDIT=dmatrix(1,1,1,Nrowmax);
	sum1=dmatrix(1,Ncov,1,Ncov);
	sum1T=dmatrix(1,Ncov,1,Ncov);
	sum2=dmatrix(1,Ncov,1,Ncov);
	sumtemp1=matrix(1, Ncov, 1, Ncov); 
	sumtemp2=matrix(1, Ncov, 1, Ncov); 
	term1=dmatrix(1,Ncov,1,Ncov);
	term2=dmatrix(1,Ncov,1,Ncov);
	varMatrixI=dmatrix(1,Nrowmax,1,Nrowmax);


    for(dll=1;dll<=Ncov;dll++)
	   for(dmm=1;dmm<=Ncov;dmm++){
			sum1[dll][dmm]=0.0;
		    sum2[dll][dmm]=0.0;
	}

	for(dii=1;dii<=Nrowtotal;dii++) 
		residualD[dii][1]=residual[dii];
 
	position=0;
	for(dii=1; dii<=Nsample; dii++){
		for(djj=1;djj<=Nrow[dii];djj++){
			position+=1;
			for(dkk=1;dkk<=Ncov;dkk++)
                designXXI[djj][dkk]=designXX[position][dkk];
             residualDI[djj][1]=residualD[position][1];
			}
		for(djj=1; djj<=Nrow[dii]; djj++)
			for(dkk=1; dkk<=Nrow[dii]; dkk++)
		 	   varMatrixI[djj][dkk]=varMatrix[dii][djj][dkk];

	   matrixtranspose(designXXI, Nrow[dii], Ncov, designXXIT);
	   matrixmultiply(designXXIT, Ncov, Nrow[dii], varMatrixI, Nrow[dii], Nrow[dii], inter3);
	   matrixmultiply(inter3, Ncov, Nrow[dii], designXXI, Nrow[dii], Ncov, term1);
  	   matrixtranspose(residualDI, Nrow[dii], 1, residualDIT);
	   matrixmultiply(residualDI, Nrow[dii], 1, residualDIT, 1, Nrow[dii], covYI);
	   matrixtranspose(inter3, Ncov, Nrow[dii], inter3T);
   	   matrixmultiply(inter3, Ncov, Nrow[dii], covYI, Nrow[dii], Nrow[dii], inter4);
	   matrixmultiply(inter4, Ncov, Nrow[dii], inter3T, Nrow[dii], Ncov, term2);
 	   for(dll=1;dll<=Ncov;dll++)
		  for(dmm=1;dmm<=Ncov;dmm++) 
			  sum1[dll][dmm]+=term1[dll][dmm];
  	   for(dll=1;dll<=Ncov;dll++)
		  for(dmm=1;dmm<=Ncov;dmm++)
				   sum2[dll][dmm]+=term2[dll][dmm];		
	}/* dii */ 

	  for(dll=1;dll<=Ncov;dll++)
		  for(dmm=1;dmm<=Ncov;dmm++)
			  sumtemp1[dll][dmm]=sum1[dll][dmm];	
	  ABSmatrixNhalf(sumtemp1, sumtemp2, Ncov, -1.0); 
	  for(dll=1;dll<=Ncov;dll++)
		  for(dmm=1;dmm<=Ncov;dmm++)
			  sum1[dll][dmm]=sumtemp2[dll][dmm];	
	  
	   matrixtranspose(sum1, Ncov, Ncov, sum1T);
       matrixmultiply(sum1, Ncov, Ncov, sum2, Ncov, Ncov, inter5);
	   matrixmultiply(inter5, Ncov, Ncov, sum1T, Ncov, Ncov, varGEEbeta);

	  free_dmatrix(covYI,1,Nrowmax,1,Nrowmax);
	  free_dmatrix(designXXI,1,Nrowmax,1,Ncov);
	  free_dmatrix(designXXIT,1,Ncov,1,Nrowmax);
	  free_dmatrix(inter3,1,Ncov,1,Nrowmax);
	  free_dmatrix(inter3T,1,Nrowmax,1,Ncov);
	  free_dmatrix(inter4,1,Ncov,1,Nrowmax);
	  free_dmatrix(inter5,1,Ncov,1,Ncov);
	  free_ivector(Nrow, 1,Nsample);
	  free_dmatrix(residualD,1,Nrowtotal,1,1);
	  free_dmatrix(residualDI,1,Nrowmax,1,1);
	  free_dmatrix(residualDIT,1,1,1,Nrowmax);
	  free_dmatrix(sum1,1,Ncov,1,Ncov);
	  free_dmatrix(sum1T,1,Ncov,1,Ncov);
	  free_dmatrix(sum2,1,Ncov,1,Ncov);
	  free_matrix(sumtemp1, 1, Ncov, 1, Ncov); 
	  free_matrix(sumtemp2, 1, Ncov, 1, Ncov); 
	  free_dmatrix(term1,1,Ncov,1,Ncov);
	  free_dmatrix(term2,1,Ncov,1,Ncov);
 	  free_dmatrix(varMatrixI,1,Nrowmax,1,Nrowmax);
  
 
}/* end */ 


void FinalCorrHongtu(float **Corr1, float **Corr2, float ***Final, int *mi, int li, int Nmax, float **ExactTime)
{

	 int dii, djj, dkk, dll, dmm;
	 float** temptMatrix=matrix(1, Nmax, 1, Nmax);
	 double* signCorr2=dvector(1, li);

	 for(dmm=1; dmm<=li; dmm++)
		 if(Corr2[dmm][dmm]>=0) signCorr2[dmm]=1.0;
		 else signCorr2[dmm]=-1.0; 

	
	 for(dmm=1; dmm<=Nsample; dmm++){
		for(dii=1; dii<=mi[dmm]; dii++)
			 for(djj=dii; djj<=mi[dmm]; djj++)
				 for(dkk=1; dkk<=li; dkk++)
					 for(dll=1; dll<=li; dll++){ 
						 if(dii!=djj) 
		 			        Final[dmm][(dii-1)*li+dkk][(djj-1)*li+dll]=Corr1[dkk][dll]*pow(fabs(Corr2[dkk][dkk])*1.0, fabs(ExactTime[dmm][dii]-ExactTime[dmm][djj])*1.0);
						 else   
						    Final[dmm][(dii-1)*li+dkk][(djj-1)*li+dll]=Corr1[dkk][dll];
                        Final[dmm][(djj-1)*li+dll][(dii-1)*li+dkk]=Final[dmm][(dii-1)*li+dkk][(djj-1)*li+dll]; 
					 }
		 for(dii=1; dii<=mi[dmm]*li; dii++){
			for(djj=1; djj<=mi[dmm]*li; djj++){ 
                temptMatrix[dii][djj]=Final[dmm][dii][djj];
				temptMatrix[djj][dii]=Final[dmm][dii][djj];
			  }
		 }

		 
		ABSmatrixNhalf(temptMatrix, Final[dmm], mi[dmm]*li, -1.0);  
	 	
//		 matrixprint(temptMatrix,mi[dmm]*li,mi[dmm]*li);

//		iv(temptMatrix, mi[dmm]*li);
//        for(dii=1; dii<=mi[dmm]*li; dii++)
//			for(djj=1; djj<=mi[dmm]*li; djj++) {
//                Final[dmm][dii][djj]=(float)temptMatrix[dii][djj];
			//	printf("\n Final=%f",Final[dmm][dii][djj]);
//			}
		 
		 
	 }/* dmm */ 

	free_matrix(temptMatrix, 1, Nmax, 1, Nmax);
	free_dvector(signCorr2, 1, li);

 }/* end */ 


double AR1_dim(double x, double** xxYmatrix, int  N2total) 
{

double fvalue=0.0;
int dii; 


for(dii=1; dii<=N2total; dii++){
//	printf("%d %lf %lf %lf \n", dii, xxYmatrix[dii][2], xxYmatrix[dii][1], pow(x, xxYmatrix[dii][1])); 
	if(fabs(xxYmatrix[dii][2])<3){ 
		if(x>=0) 
         fvalue+=(xxYmatrix[dii][2]-pow(x, xxYmatrix[dii][1]))*(xxYmatrix[dii][2]-pow(x, xxYmatrix[dii][1])); 
		else 
         fvalue+=(xxYmatrix[dii][2]+pow(fabs(x), xxYmatrix[dii][1]))*(xxYmatrix[dii][2]+pow(fabs(x), xxYmatrix[dii][1])); 
	}
}


return(fvalue); 

}/* f1dim */  




double goldenHT(double ax, double bx, double cx, double (*f)(double, double**, int),  double tol,
	double *xmin, double **xxY, int N2total)
{
   double f1,f2,x0,x1,x2,x3;
   double R=0.61803399; 
   double C=(1.0-R);
   int iter=1; 

	x0=ax;
	x3=cx;
	if (fabs(cx-bx) > fabs(bx-ax)) {
		x1=bx;
		x2=bx+C*(cx-bx);
	} else {
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=(*f)(x1, xxY, N2total);
	f2=(*f)(x2, xxY, N2total);
	while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))&&iter<=100) {
		if (f2 < f1) {
			SHFT3(x0,x1,x2,R*x1+C*x3)
			SHFT2(f1,f2,(*f)(x2, xxY, N2total))
		} else {
			SHFT3(x3,x2,x1,R*x2+C*x0)
			SHFT2(f2,f1,(*f)(x1, xxY, N2total))
		}
		iter++; 
	}
	if (f1 < f2) {
		*xmin=x1;
		return f1;
	} else {
		*xmin=x2;
		return f2;
	}
}



double EstimateAR_rho(double** timeDresidD, int* indxMi, int djj, int N2total) 
{
   int dii, dmm, dnn, dkk, dll; 
   double** xxYmatrix=dmatrix(1, N2total, 1, 2); 
   double xx,xmin,fx,fb,fa,bx,ax, cx;

   for(dii=1; dii<=N2total; dii++){
      xxYmatrix[dii][1]=timeDresidD[dii][1];     
      xxYmatrix[dii][2]=timeDresidD[dii][djj+1];   
   } 
   ax=0.0;
   bx=0.25;
   cx=0.99; 
//   mnbrakHT(&ax, &cx, &bx, &fa, &fx, &fb,AR1_dim, xxYmatrix, N2total);
   fa=1.0e-3;

   fx=goldenHT(ax,  bx,  cx, AR1_dim,  fa, &fb,  xxYmatrix,  N2total); 
   free_dmatrix(xxYmatrix, 1, N2total, 1, 2); 
   return(fb);
 
}/* end */ 




void	CalVarWeightedGEECorrectedForVar(int current,double **outputYYO,double **designXX,int li,int *mi,float **beta,int q, double **varGEE,double *weight, int *labelid, int noneighs, float ****VarMatrix)
{

	double  **designXXI,**designXXIT,**inter1,**inter2,**inter3,**inter4,**sum1,**sum2,**sum2T,**newsum2,**sum3,**sum4,**outputYYOD,**residualD,**residualDI,**residualDIT,**temp1,**temp1T,**temp2,**temp3,**temp4,**VarMatrixD;
    int dii,djj,dkk,dll,dmm, Ntotal,position,Npar,mimax,Nrowmax,*Nrow,Nrowtotal,rankXX,NoP0,NcovVar; 
	float **HONGterm1, **HONGterm2; 


	mimax=0;
	Nrowmax=0;

	//printf("\n Nsample=%d",Nsample);
	
	Nrow=ivector(1,Nsample);
	Nrowtotal=0;


	for (dii=1;dii<=Nsample;dii++){
		Nrow[dii]=mi[dii]*li;
	//	Nrow[dii]=3;
	
	if (mi[dii]>mimax) mimax=mi[dii];
	//	mimax=1;

	    Nrowtotal+=Nrow[dii];
	}



     Nrowmax=mimax*li;
	 NcovVar=li*q;
	// Ncov=2;

//	 printf("\n Nrowmax=%d,Ncov=%d",Nrowmax,Ncov);

	designXXI=dmatrix(1,Nrowmax,1,Ncov);
	designXXIT=dmatrix(1,Ncov,1,Nrowmax);
	HONGterm1=matrix(1, Ncov, 1, Ncov); 
	HONGterm2=matrix(1, Ncov, 1, Ncov); 
	inter1=dmatrix(1,Ncov,1,Nrowmax);
    inter2=dmatrix(1,Ncov,1,Ncov);
	inter3=dmatrix(1,Ncov,1,Nrowmax);
	inter4=dmatrix(1,Ncov,1,Ncov);

	sum1=dmatrix(1,Ncov,1,Ncov);
	sum2=dmatrix(1,Nrowmax,1,Nrowmax);
	newsum2=dmatrix(1,Nrowmax,1,1);
	sum2T=dmatrix(1,1,1,Nrowmax);
	sum3=dmatrix(1,Ncov,1,Ncov);
	sum4=dmatrix(1,Ncov,1,Ncov);

	outputYYOD=dmatrix(1,Nrowmax,1,1);
	residualD=dmatrix(1,Nrowtotal,1,1);
	residualDI=dmatrix(1,Nrowmax,1,1);
	residualDIT=dmatrix(1,1,1,Nrowmax);
	
		temp1=dmatrix(1,Nrowmax,1,1);
		temp1T=dmatrix(1,1,1,Nrowmax);
		temp2=dmatrix(1,Nrowmax,1,Nrowmax);
		temp3=dmatrix(1,Ncov,1,Ncov);
		temp4=dmatrix(1,Nrowmax,1,Nrowmax);
	double **temp1a=dmatrix(1,Ncov,1,Nrowmax);
	double **temp2a=dmatrix(1,Ncov,1,Ncov);

	VarMatrixD=dmatrix(1,Nrowmax,1,Nrowmax);

	double newsum1=0;
		

	for (dll=1;dll<=Nrowmax;dll++){
		newsum2[dll][1]=0.0;
	for (dmm=1;dmm<=Nrowmax;dmm++){
			temp2[dll][dmm]=0.0;
			//temp3[dll][dmm]=0.0;
			sum2[dll][dmm]=0.0;
			}
	}



for (dll=1;dll<=Ncov;dll++)
	for (dmm=1;dmm<=Ncov;dmm++){
			sum4[dll][dmm]=0.0;
			sum3[dll][dmm]=0.0;
			temp3[dll][dmm]=0.0;
			}


	for (dll=1;dll<=Nrowmax;dll++){
          temp1[dll][1]=0.0;
		  temp1T[1][dll]=0.0;
		  }


int dnn,ii,jj;
int pos=0;

	for (dii=1; dii<=Nsample; dii++){
		newsum1=0;
		
		for (dll=1;dll<=Nrowmax;dll++){
			newsum2[dll][1]=0.0;
			for (dmm=1;dmm<=Nrowmax;dmm++){
				sum2[dll][dmm]=0.0;
			}
		}

		for(dll=1;dll<=Ncov;dll++)
			for(dmm=1;dmm<=Ncov;dmm++)
				sum1[dll][dmm]=0;


		for (djj=1;djj<=Nrow[dii];djj++){
              pos+=1;
			for (dkk=1;dkk<=Ncov;dkk++)
                 designXXI[djj][dkk]=designXX[pos][dkk];
			

		}


		matrixtranspose(designXXI,Nrow[dii],Ncov,designXXIT);


			for (dnn=1;dnn<=noneighs;dnn++){  
				for (djj=1;djj<=Nrow[dii];djj++){			  
					 outputYYOD[djj][1]=outputYYO[labelid[dnn]][pos-Nrow[dii]+djj];
					 residualDI[djj][1]=0.0;            
			        for (dkk=1;dkk<=Ncov;dkk++)
				       residualDI[djj][1]+=designXXI[djj][dkk]*beta[labelid[dnn]][dkk];
  			        residualDI[djj][1]=outputYYOD[djj][1]-residualDI[djj][1];
		   }//djj

  		   for (djj=1;djj<=Nrow[dii];djj++)
				for(dkk=1;dkk<=Nrow[dii];dkk++){
						VarMatrixD[djj][dkk]=VarMatrix[labelid[dnn]][dii][djj][dkk];
					}

		    matrixmultiply(VarMatrixD,Nrow[dii],Nrow[dii],residualDI,Nrow[dii],1,temp1);
			matrixmultiply(designXXIT,Ncov,Nrow[dii],VarMatrixD,Nrow[dii],Nrow[dii],temp1a);
            matrixmultiply(temp1a,Ncov,Nrow[dii],designXXI,Nrow[dii],Ncov,temp2a);


			    for (ii=1;ii<=Ncov;ii++)
					for(jj=1;jj<=Ncov;jj++){
						sum1[ii][jj]+=weight[dnn]*temp2a[ii][jj];
					}

				for (ii=1;ii<=Nrow[dii];ii++){
					newsum2[ii][1]+=weight[dnn]*temp1[ii][1];

					}
					
	}//dnn,Noneighs



    matrixtranspose(newsum2,Nrow[dii],1,sum2T);
  	matrixmultiply(newsum2,Nrow[dii],1,sum2T,1,Nrow[dii],sum2);
    matrixmultiply(designXXIT,Ncov,Nrow[dii],sum2,Nrow[dii],Nrow[dii],inter3);
	matrixmultiply(inter3,Ncov,Nrow[dii],designXXI,Nrow[dii],Ncov,inter4);

 
	for (dll=1;dll<=Ncov;dll++)
		for (dmm=1;dmm<=Ncov;dmm++)
			sum3[dll][dmm]+=sum1[dll][dmm];
 
		
	for (dll=1;dll<=Ncov;dll++)
		for (dmm=1;dmm<=Ncov;dmm++)
			sum4[dll][dmm]+=inter4[dll][dmm];
 

	}//dii,Nsample


 
//		iv(sum3,Ncov);
	
	for (djj=1;djj<=Ncov;djj++)
		for (dkk=1;dkk<=Ncov;dkk++)
			HONGterm1[djj][dkk]=sum3[djj][dkk]; 
	
	//	iv(term1,Ncov);
	
	ABSmatrixNhalf(HONGterm1, HONGterm2, Ncov, -1.0); 
	
	for (djj=1;djj<=Ncov;djj++)
		for (dkk=1;dkk<=Ncov;dkk++)
			sum3[djj][dkk]=HONGterm2[djj][dkk]; 
	

		matrixmultiply(sum3,Ncov,Ncov,sum4,Ncov,Ncov,temp3);
		matrixmultiply(temp3,Ncov,Ncov,sum3,Ncov,Ncov,varGEE);

	

//	matrixprint(varGEE,Ncov,Ncov);
	

	free_ivector(Nrow,1,Nsample);	   
 	free_dmatrix(designXXI,1,Nrowmax,1,Ncov);
 	free_dmatrix(designXXIT,1,Ncov,1,Nrowmax);
	free_matrix(HONGterm1, 1, Ncov, 1, Ncov); 
	free_matrix(HONGterm2, 1, Ncov, 1, Ncov); 
	free_dmatrix(inter1,1,Ncov,1,Nrowmax);
    free_dmatrix(inter2,1,Ncov,1,Ncov);
	free_dmatrix(inter3,1,Ncov,1,Nrowmax);
	free_dmatrix(inter4,1,Ncov,1,Ncov);
	free_dmatrix(sum1,1,Ncov,1,Ncov);
	free_dmatrix(sum2,1,Nrowmax,1,Nrowmax);
	free_dmatrix(newsum2,1,Nrowmax,1,1);
	free_dmatrix(sum2T,1,1,1,Nrowmax);
	free_dmatrix(sum3,1,Ncov,1,Ncov);
	free_dmatrix(sum4,1,Ncov,1,Ncov);
	free_dmatrix(outputYYOD,1,Nrowmax,1,1);
	free_dmatrix(residualD,1,Nrowtotal,1,1);
	free_dmatrix(residualDI,1,Nrowmax,1,1);
	free_dmatrix(residualDIT,1,1,1,Nrowmax);
	free_dmatrix(temp1,1,Nrowmax,1,1);
	free_dmatrix(temp1T,1,1,1,Nrowmax);
	free_dmatrix(temp2,1,Nrowmax,1,Nrowmax);
	free_dmatrix(temp3,1,Ncov,1,Ncov);
	free_dmatrix(temp4,1,Nrowmax,1,Nrowmax);
	
	free_dmatrix(temp1a,1,Ncov,1,Nrowmax);
	free_dmatrix(temp2a,1,Ncov,1,Ncov);

	free_dmatrix(VarMatrixD,1,Nrowmax,1,Nrowmax);


}

void CalMeanParaHongtu(double **designXX, int p, float ***varMatrix,   double *residual, double *beta, int *mi, int li,double *outputYYO, double *oldbeta)
{
		
	double  *betanew,**inter1,**inter2,**inter3,**sum1,**sum2,*sum3,**term1,**term2,**term3; 
	double   **designXXI,**designXXIT,  **residualD, **residualDI;
    int dii, djj, dkk, dll, dmm, Ntotal, position, Npar, mimax, Nrowmax, *Nrow, Nrowtotal, rankXX, NoP0; 
	float **sumtemp1, **sumtemp2; 
	
	double **varMatrix2;
	mimax=0;
	Nrowmax=0;
	
	Nrow=ivector(1,Nsample);
	Nrowtotal=0;
	for (dii=1;dii<=Nsample;dii++){
		Nrow[dii]=mi[dii]*li;
 	    if(mi[dii]>mimax) mimax=mi[dii];
	    Nrowtotal+=Nrow[dii];
	}
 
	Nrowmax=mimax*li;


	betanew=dvector(1,Ncov);
	designXXI=dmatrix(1,Nrowmax,1,Ncov);
	designXXIT=dmatrix(1,Ncov,1,Nrowmax);
	inter1=dmatrix(1,Nrowmax,1,Nrowmax);
    inter2=dmatrix(1,Nrowmax,1,Nrowmax);
	inter3=dmatrix(1,Ncov,1,Nrowmax);
	residualD=dmatrix(1,Nrowtotal,1,1);
	residualDI=dmatrix(1,Nrowmax,1,1);
	sum1=dmatrix(1,Ncov,1,Ncov);
	sum2=dmatrix(1,Ncov,1,1);
	sum3=dvector(1,Nrowtotal);
	sumtemp1=matrix(1, Ncov, 1, Ncov); 
	sumtemp2=matrix(1, Ncov, 1, Ncov); 
	term1=dmatrix(1,Ncov,1,Ncov);
	term2=dmatrix(1,Ncov,1,1);
	term3=dmatrix(1,Ncov,1,1);
	varMatrix2=dmatrix(1, Nrowmax, 1, Nrowmax); 
 		
	for(dii=1;dii<=Nrowtotal;dii++)
		residualD[dii][1]=residual[dii];
 	for(dll=1;dll<=Ncov;dll++)
		for (dmm=1;dmm<=Ncov;dmm++)
			sum1[dll][dmm]=0.0;
	for(dll=1;dll<=Ncov;dll++)
		sum2[dll][1]=0.0;
    position=0;
	for(dii=1; dii<=Nsample; dii++){
		for(djj=1; djj<=Nrow[dii]; djj++)
			for(dkk=1; dkk<=Ncov; dkk++)
                designXXI[djj][dkk]=0.0;
   	    for(djj=1; djj<=Nrow[dii]; djj++){
			position+=1;
			for(dkk=1; dkk<=Ncov; dkk++) 
                designXXI[djj][dkk]=designXX[position][dkk];
            residualDI[djj][1]=residualD[position][1];
			}
		for(djj=1; djj<=Nrow[dii]; djj++)
			for(dkk=1; dkk<=Nrow[dii]; dkk++)
                varMatrix2[djj][dkk]=varMatrix[dii][djj][dkk]*1.0;  
		matrixtranspose(designXXI, Nrow[dii], Ncov, designXXIT);
	    matrixmultiply(designXXIT, Ncov, Nrow[dii], varMatrix2, Nrow[dii], Nrow[dii], inter3);
	    matrixmultiply(inter3, Ncov, Nrow[dii], designXXI, Nrow[dii], Ncov, term1);
	    matrixmultiply(inter3, Ncov, Nrow[dii], residualDI, Nrow[dii], 1, term2);
 	    for(dll=1;dll<=Ncov;dll++)
		   for(dmm=1;dmm<=Ncov;dmm++)
			  sum1[dll][dmm]+=term1[dll][dmm];
	    for(dll=1;dll<=Ncov;dll++)
			 sum2[dll][1]+=term2[dll][1];
		
	}/* for */ 

   for(dii=1; dii<=Ncov; dii++)
		for(djj=1; djj<=Ncov; djj++)	
			sumtemp1[dii][djj]=sum1[dii][djj];
	
   ABSmatrixNhalf(sumtemp1, sumtemp2, Ncov, -1.0); 
   for(dii=1; dii<=Ncov; dii++)
	   for(djj=1; djj<=Ncov; djj++)	
		  sum1[dii][djj]=sumtemp2[dii][djj];	
	
   matrixmultiply(sum1, Ncov, Ncov, sum2, Ncov, 1, term3);
   for(dii=1;dii<=Ncov;dii++) 
     betanew[dii]=beta[dii]+term3[dii][1];

   for(dii=1;dii<=Ncov;dii++)
	  oldbeta[dii]=beta[dii];
   for(dii=1;dii<=Ncov;dii++){
   	  beta[dii]=betanew[dii];
//	  printf("\nbeta[dii]=%f",beta[dii]);
	  }

   for(dii=1;dii<=Nrowtotal;dii++)
	   sum3[dii]=0.0;
   for(dii=1;dii<=Nrowtotal;dii++)
	   for(djj=1;djj<=Ncov;djj++)
		   sum3[dii]+=designXX[dii][djj]*beta[djj];
   for(dii=1;dii<=Nrowtotal;dii++)
  	  residual[dii]=outputYYO[dii]-sum3[dii];



	free_dvector(betanew,1,Ncov);
	free_dmatrix(designXXI,1,Nrowmax,1,Ncov);
	free_dmatrix(designXXIT,1,Ncov,1,Nrowmax);
	free_dmatrix(inter1,1,Nrowmax,1,Nrowmax);
    free_dmatrix(inter2,1,Nrowmax,1,Nrowmax);
	free_dmatrix(inter3,1,Ncov,1,Nrowmax);
	free_ivector(Nrow,1,Nsample);
	free_dmatrix(residualD,1,Nrowtotal,1,1);
	free_dmatrix(residualDI,1,Nrowmax,1,1);
	free_dmatrix(sum1,1,Ncov,1,Ncov);
	free_dmatrix(sum2,1,Ncov,1,1);
	free_dvector(sum3,1,Nrowtotal);
	free_matrix(sumtemp1, 1, Ncov, 1, Ncov); 
	free_matrix(sumtemp2, 1, Ncov, 1, Ncov); 
	free_dmatrix(term1,1,Ncov,1,Ncov);
	free_dmatrix(term2,1,Ncov,1,1);
	free_dmatrix(term3,1,Ncov,1,1);
	free_dmatrix(varMatrix2, 1, Nrowmax, 1, Nrowmax); 


}/* end */ 



void GEEestimatesHongtu(double *beta,   double *residual, double **designXX, double **TXX, double *outputYYO, int li, int *mi, int p, int q, float ***varMatrix, float **ExactTime, int* indxMI, int* indxMI1)
{
	double  diff, *oldbeta, *oldxi;
	float  **Corr1,**Corr2;
	int dii, dkk, djj, dmm, rankXX,  NoP0,  obsforonemax,mimax,  *Nrow,  Nrowtotal;

    Nrow=ivector(1,Nsample);
	Nrowtotal=0;
    mimax=0;
 
	for (dii=1;dii<=Nsample;dii++){
		Nrow[dii]=mi[dii]*li;
		if (mi[dii]>mimax)
			mimax=mi[dii];
 		Nrowtotal+=Nrow[dii];
	}

	obsforonemax=mimax*li;
	Corr1=matrix(1,li,1,li);
	Corr2=matrix(1,li,1,li);
 	oldbeta=dvector(1,Ncov);

//	printf("\n obsforonemax=%d",obsforonemax);
 
    linearHongtu(beta,  residual, designXX,  TXX, outputYYO);   

	//	for(dii=1; dii<=Ncov; dii++){
	  // printf("\n beta=%f",beta[dii]);
	 //  }


	for (dii=1;dii<=Nsample;dii++)
		for (djj=1;djj<=obsforonemax;djj++)
			for (dkk=1;dkk<=obsforonemax;dkk++)
				varMatrix[dii][djj][dkk]=0.0;
	diff=0.0;
    dmm=1;

	do{
      AR1TimeHongtu(residual, mi, li, Corr1, Corr2, ExactTime, indxMI,  indxMI1);   
 	  FinalCorrHongtu(Corr1, Corr2, varMatrix, mi, li, mimax*li, ExactTime);//need to check, which way is correct to do the product
  	  CalMeanParaHongtu(designXX, p, varMatrix,  residual, beta, mi, li, outputYYO, oldbeta);
	  dmm++;
	}while(dmm<=3);

 

	for(dii=Ncov+1; dii<=Ncov+DimSPD; dii++){
       beta[dii]=Corr1[dii-Ncov][dii-Ncov]; 
	 //  printf("\n beta=%f",beta[dii]);
	   }

    for(dii=Ncov+DimSPD+1; dii<=Ncov+2*DimSPD; dii++) {
       beta[dii]=Corr2[dii-Ncov-DimSPD][dii-Ncov-DimSPD]; 
     //  printf("\n beta=%f",beta[dii]);
	   }

	
 
    free_ivector(Nrow,1,Nsample);
	free_matrix(Corr1, 1, li, 1, li);
	free_matrix(Corr2, 1,li,1,li);
	free_dvector(oldbeta,1,Ncov);
}



float waldtest(float *BrainTheta, double **varGEE,float **RR,float *rr00,int noRow)

{

	int dii,djj,dkk,dll;
	float* temp1=vector(1, noRow);
	float** temp2=matrix(1, noRow, 1,Ncov);
	float** temp3=matrix(1, noRow,1, noRow);
	float* temp4=vector(1, noRow);
    float waldtest;

 
    for(dii=1;dii<=noRow;dii++)
		 temp4[dii]=0.0;
	for(dii=1;dii<=noRow;dii++){
		temp1[dii]=0.0;
		for(djj=1;djj<=Ncov;djj++) 
			temp1[dii]+=RR[dii-1][djj-1]*BrainTheta[djj];
		temp1[dii]-=rr00[dii-1];
	} 
    for(dii=1; dii<=noRow; dii++)
		for(djj=1; djj<=Ncov; djj++){
			temp2[dii][djj]=0.0; 
			for(dkk=1; dkk<=Ncov; dkk++)
			   temp2[dii][djj]+=RR[dii-1][dkk-1]*varGEE[djj][dkk];
		} 
   for(dii=1; dii<=noRow; dii++)
		for(djj=1; djj<=noRow; djj++){
            temp3[dii][djj]=0.0;
			for(dkk=1; dkk<=Ncov; dkk++)
				temp3[dii][djj]+=temp2[dii][dkk]*RR[djj-1][dkk-1];
		} 
	ABSmatrixNhalf(temp3, temp2, noRow, -1.0);  
	waldtest=0.0;
     for (dii=1; dii<=noRow; dii++)
		 for(djj=1; djj<=noRow; djj++)
		     waldtest+=temp1[djj]*temp2[djj][dii]*temp1[dii];
	
	free_vector(temp1,1,noRow);
	free_matrix(temp2,1,noRow,1,Ncov);
	free_matrix(temp3,1,noRow,1,noRow);
	free_vector(temp4,1,noRow);
	return(waldtest);

}



double gammln(double xx)
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}/*  gammln */ 




double  betacf(double a,double b,double x)
{

    double aa,c,d,del,h,qab,qam,qap;
    double epSilon=3.0e-7; 
	int m, m2, NitMax=100; 
    double FPMIN=1.0e-30;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;  
    d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=NitMax;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;  
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;  
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < epSilon) break;  
//		printf("%14lf ", h); 
		}
//    if (m > NitMax) 
//        nrerror("a or b too big, or MAXIT too small in betacf");
    return h;

}/*betacf*/ 
  



double betai(double a, double b, double x)
{
	double bt;
 
	 
	if (x < 0.0 || x > 1.0) 
		nrerror("Bad x in routine BETAI");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x);
 	if (x < (a+1.0)/(a+b+2.0))
		return exp(bt)*betacf(a,b,x)/a;
	else
		return 1.0-exp(bt)*betacf(b,a,1.0-x)/b;
}/* betai */ 


void CalMeanParaWeightGEECorrectedForVar (double **designXX, float *beta,int *mi, int li,double **outputYYO,double *weight,int *labelid,int noneighs,float ****VarMatrix)

{
	
	
	double  **term1,**term2,**term3,**designXXI,**designXXIT,**VarMatrixD,**temp1,**temp2,**temp3;
    int dii,djj,dkk,dll,dmm, Ntotal,position,Npar,mimax,Nrowmax,*Nrow,Nrowtotal,rankXX,NoP0; 
	double **outputYYOD;
	float **HONGterm1, **HONGterm2; 


	mimax=0;
	Nrowmax=0;
	
	Nrow=ivector(1,Nsample);
	Nrowtotal=0;

	for (dii=1;dii<=Nsample;dii++){
		Nrow[dii]=mi[dii]*li;
	
	if (mi[dii]>mimax) mimax=mi[dii];
	    Nrowtotal+=Nrow[dii];
	//	printf("\n Nsample=%d, Ncov=%d, Nrow[dii]=%d  mi[dii]=%d",Nsample,Ncov, Nrow[dii],mi[dii]);
	}

     Nrowmax=mimax*li;

//	printf("\n Nrowmax=%d, Nrowtotal=%d",Nrowmax,Nrowtotal);

	
	designXXI=dmatrix(1,Nrowmax,1,Ncov);
	designXXIT=dmatrix(1,Ncov,1,Nrowmax);
	HONGterm1=matrix(1, Ncov, 1, Ncov); 
	HONGterm2=matrix(1, Ncov, 1, Ncov); 
	temp1=dmatrix(1,Ncov,1,Nrowmax);
    temp2=dmatrix(1,Ncov,1,Ncov);
	temp3=dmatrix(1,Ncov,1,1);
	term1=dmatrix(1,Ncov,1,Ncov);
	term2=dmatrix(1,Ncov,1,1);
	term3=dmatrix(1,Ncov,1,1);
	outputYYOD=dmatrix(1,Nrowmax,1,1);
	VarMatrixD=dmatrix(1,Nrowmax,1,Nrowmax);

  //   printf("\n where is wrong 2?");


	for (dll=1;dll<=Ncov;dll++)
		for (dmm=1;dmm<=Nrowmax;dmm++)
			temp1[dll][dmm]=0.0;


	for (dll=1;dll<=Ncov;dll++)
		for (dmm=1;dmm<=Ncov;dmm++){
		    temp2[dll][dmm]=0.0;
			term1[dll][dmm]=0.0;
		}


	for (dll=1;dll<=Ncov;dll++){
		temp3[dll][1]=0.0;
		term2[dll][1]=0.0;
		term3[dll][1]=0.0;
		
	}


 
	int dnn,start=0;
	double sumweight=0.0,temp;

	for (dnn=1;dnn<=noneighs;dnn++){
           start=0;
		for (dii=1;dii<=Nsample;dii++){

            
			
			for (djj=1;djj<=Nrow[dii];djj++){
				outputYYOD[djj][1]=outputYYO[labelid[dnn]][start+djj];
			//	printf("\n outputYYOD[djj][1]=%f",outputYYOD[djj][1]);

				for (dkk=1;dkk<=Ncov;dkk++){
				//	printf("\nstart+djj=%d,designXX[start+djj][dkk]=%f",start+djj,designXX[start+djj][dkk]);
					designXXI[djj][dkk]=designXX[start+djj][dkk];
			//		printf("\n dnn=%d,designXXI[djj][dkk]=%f",dnn,designXXI[djj][dkk]);
			}
			}

				start+=Nrow[dii];

				matrixtranspose(designXXI,Nrow[dii],Ncov,designXXIT);

				for (djj=1;djj<=Nrow[dii];djj++)
					for (dkk=1;dkk<=Nrow[dii];dkk++)
						VarMatrixD[djj][dkk]=VarMatrix[labelid[dnn]][dii][djj][dkk];


                   //  matrixprint(VarMatrixD,Nrow[dii],Nrow[dii]);
				//	iv(VarMatrixD,Nrow[dii]);// Yimei Sept 7th.2009

					

			matrixmultiply(designXXIT,Ncov,Nrow[dii],VarMatrixD,Nrow[dii],Nrow[dii],temp1);
		//	matrixprint(temp1,Ncov,Nrow[dii]);
			matrixmultiply(temp1,Ncov,Nrow[dii],designXXI,Nrow[dii],Ncov,temp2);
		//	matrixmultiply(designXXIT,Ncov,Nrow[dii],designXXI,Nrow[dii],Ncov,temp2);

			for (djj=1;djj<=Ncov;djj++)
				for (dkk=1;dkk<=Ncov;dkk++)
					term1[djj][dkk]+=weight[dnn]*temp2[djj][dkk];


			matrixmultiply(temp1,Ncov,Nrow[dii],outputYYOD,Nrow[dii],1,temp3);
		//	matrixmultiply(designXXIT,Ncov,Nrow[dii],outputYYOD,Nrow[dii],1,temp3);


            for (djj=1;djj<=Ncov;djj++)
				term2[djj][1]+=weight[dnn]*temp3[djj][1];

		}//dii


		
	}//dnn

	for (djj=1;djj<=Ncov;djj++)
		for (dkk=1;dkk<=Ncov;dkk++)
			HONGterm1[djj][dkk]=term1[djj][dkk]; 
			
//	iv(term1,Ncov);
	
	ABSmatrixNhalf(HONGterm1, HONGterm2, Ncov, -1.0); 

	for (djj=1;djj<=Ncov;djj++)
		for (dkk=1;dkk<=Ncov;dkk++)
			term1[djj][dkk]=HONGterm2[djj][dkk]; 
	
	matrixmultiply(term1,Ncov,Ncov,term2,Ncov,1,term3);



		for (dii=1;dii<=Ncov;dii++){
			beta[dii]=term3[dii][1];
		//	printf("\n%f",beta[dii]);
		}




	free_ivector(Nrow,1,Nsample);
	free_dmatrix(designXXI,1,Nrowmax,1,Ncov);
	free_dmatrix(designXXIT,1,Ncov,1,Nrowmax);
	free_matrix(HONGterm1, 1, Ncov, 1, Ncov); 
	free_matrix(HONGterm2, 1, Ncov, 1, Ncov); 
	free_dmatrix(temp1,1,Ncov,1,Nrowmax);
    free_dmatrix(temp2,1,Ncov,1,Ncov);
	free_dmatrix(temp3,1,Ncov,1,1);
	free_dmatrix(term1,1,Ncov,1,Ncov);
	free_dmatrix(term2,1,Ncov,1,1);
	free_dmatrix(term3,1,Ncov,1,1);
	free_dmatrix(outputYYOD,1,Nrowmax,1,1);
	free_dmatrix(VarMatrixD,1,Nrowmax,1,Nrowmax);



}




int dmin( int a, int b){
	if (a>b)
		return(b);
	else return(a);
}


int dmax( int a, int b){
	if (a>b)
		return(a);
	else return(b);
}



void choldcYS (float **a, int n, float p[])
{

   int i,j,k;
   float sum;

  for(i=1; i<=n; i++){
	 for(j=i; j<=n;j++){
	    for(sum=a[i][j], k=i-1; k>=1; k--)  sum-=a[i][k]*a[j][k];
	       if (i==j){
	         if(sum<=0.0)
		        printf(" choldc failed");
	         p[i]=sqrt(sum);
	  } else a[j][i]=sum/p[i];
	 }
   }
 }


void cholslYS(float **a,int n,float p[],float b[],float x[])
  {
     int i, k;
     float sum;

     for (i=1; i<=n;i++){
	  for(sum=b[i], k=i-1; k>=1;k--)  sum-=a[i][k]*x[k];
	     x[i]=sum/p[i];
      }
      for (i=n; i>=1; i--){
	     for (sum=x[i], k=i+1; k<=n; k++)  sum-=a[k][i]*x[k];
	    x[i]=sum/p[i];
	  }

}


void ivYS(float **W, int Nlow, int Nup)
{
       int  i, k;
       float *XX, *PP, *COL, **WW, **WW1;
	   int NP=Nup-Nlow+1;

	   COL=vector(1, NP+1);
	   PP=vector(1, NP+1);
	   WW=matrix(1, NP+1, 1, NP+1);
	   WW1=matrix(1, NP+1, 1, NP+1);
	   XX=vector(1, NP+1);

	   for(i=Nlow; i<=Nup; i++)
		   for(k=Nlow; k<=Nup; k++)
			   WW1[i-Nlow+1][k-Nlow+1]=W[i][k];

       choldcYS(WW1, NP, PP);
       for(k=1;k<=NP;k++){
	      for(i=1;i<=NP;i++)
	         COL[i]=0.0;
	      COL[k]=1.0;
	      cholslYS(WW1, NP, PP, COL, XX);
	      for(i=1; i<=NP; i++)
			  WW[i][k]=XX[i];
	  }
	  for(i=Nlow; i<=Nup; i++)
		 for(k=Nlow; k<=Nup; k++)
		     W[i][k]=WW[i-Nlow+1][k-Nlow+1];

	   free_vector(COL, 1, NP+1);
	   free_vector(PP, 1, NP+1);
	   free_vector(XX, 1, NP+1);
	   free_matrix(WW, 1, NP+1, 1, NP+1);
       free_matrix(WW1, 1, NP+1, 1, NP+1);

 }/* end of inver.c */