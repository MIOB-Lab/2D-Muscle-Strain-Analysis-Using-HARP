#include <mex.h>

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>

#define pi		3.1415926

#define WRAP(x)		((x) - 2*pi*floor(((x)+pi)/(2*pi)))
#define SIGN(x)		((fabs(x)!= 0? 1: 0) * (x>=0? 1:-1))
#define swap(t, x, y)       do {t z=x; x=y; y=z;} while(0)

// double linked List
typedef struct SSDL_node_def {
    int x;
    int y;
	int region_label;
    double accumulatedCost;
//	double val;
//  double mx;
//  double my;
	// means the closest point that 
	// make the total cost of this point smallest
	int closest_x;	
	int closest_y;
    struct SSDL_node_def *next;
	struct SSDL_node_def *prev;
} SSDL_node;

typedef struct RegionGrowNode {
    int x;
    int y;
    struct RegionGrowNode *next;
} RG_node;


double BilinearInterpolation(double **mat, int XO, int YO, double x, double y);

double BilinearInterpolationWrap(double **mat, int sx, int sy, double x, double y);
double BilinearInterpolationWrap2(double **mat, int sx, int sy, double x, double y, double *grad);

double **pgDmatrix(int nrl,int nrh,int ncl,int nch);
void pgFreeDmatrix(double **m,int nrl,int nrh,int ncl,int nch);
int **pgImatrix(int nrl,int nrh,int ncl,int nch);
void pgFreeImatrix(int **m,int nrl,int nrh,int ncl,int nch);
unsigned long long **pgULmatrix(int nrl,int nrh,int ncl,int nch);
void pgFreeULmatrix(unsigned long long **m,int nrl,int nrh,int ncl,int nch);

int pgWriteDouble(double *ptr, int size, const char *filename);

int pgWriteDouble(double *ptr, int size, const char *filename);




// some utilities
void QKSort2(int n, double* arr, int *brr);
void QKSort3(int n, double* arr, int *brr, int *crr);

SSDL_node* InsertSSDL(SSDL_node *start, SSDL_node *newNode);
int SSDLLength(SSDL_node *start);


void TrackingRefinement_FastMarching(double **phase1_t0, double **phase1_t1, double **phase2_t0, double **phase2_t1,
						double **weight1, double **weight2,
                        int nx, int ny,
                        double *startPt,
						double *trackedPt,
                        double **motion_x,
                        double **motion_y, 
                        double **trackIndex,
						double **savedWeight, double **savedWeight1,double **savedWeight2
                        );


int TrackPoint(double **phase1_t1, double **phase2_t1, 
				  int nx, int ny,
				  double* point2D, 
				  double* currentPhase,
				  double* trackedPt,
                  double* motion
				  );

double **IdentifyGoodRegion(double **motion_x, double **motion_y, 
							double **trackIndex, double *startPt, int nx, int ny);


