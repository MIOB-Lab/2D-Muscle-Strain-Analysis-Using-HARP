
#include "HarpTrackingRefinement_FastMarching.h"


/* 
	to call:
	[return] = HarpTrackingRefinement(seed, phase1_t0, phase1_t1, phase2_t0, phase2_t1, weight1, weight2) 
	or:SSDL_node
	[return] = HarpTrackingRefinement(seed, phase1_t0, phase1_t1, phase2_t0, phase2_t1, weight1, weight2, seedInSecondFrame)

*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
	
	double startPoint[2];
	double **phase1_t0, **phase1_t1, **phase2_t0, **phase2_t1;
	double* pd;
    
    const long long unsigned int *dim_array;
    
	int nx, ny, i, j;
    int ret;
    
    double ptPhase[2];
	double ptMotion[2];
    double trackedPt[2];
    double motion[2];
    double **motion_x, **motion_y;
    double **trackIndex;
    double **rgIndex;
	
	double **weight1, **weight2;

	double maxWeight1, maxWeight2;
	double **savedWeight, **savedWeight1, **savedWeight2; // for debug purpose
	
	if (nrhs < 7)
		mexErrMsgTxt("At least 7 input variables are required!\n Usage: \n"
				"[return] = HarpTrackingRefinement_FastMarching(seed, phase1_t0, phase1_t1, phase2_t0, phase2_t1, weight1, weight2)\n"
				"or: \n"
				"[return] = HarpTrackingRefinement_FastMarching(seed, phase1_t0, phase1_t1, phase2_t0, phase2_t1, weight1, weight2, seedInSecondFrame)");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) 
		|| !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[5]) || !mxIsDouble(prhs[6]))
		mexErrMsgTxt("Input array must be of type double.");
		
	if (nrhs > 7 )
	{
		if (!mxIsDouble(prhs[7]))
			mexErrMsgTxt("Input array must be of type double.");
	}


	pd = (double*)mxGetPr(prhs[0]);
	startPoint[0] = pd[0]-1;
	startPoint[1] = pd[1]-1;

	
	dim_array = mxGetDimensions(prhs[1]);
	nx = dim_array[0];
	ny = dim_array[1];
	dim_array = mxGetDimensions(prhs[2]);
	if (nx != dim_array[0] || ny != dim_array[1]) {
		mexErrMsgTxt("The size of the phase images are not equal.");
	}
	dim_array = mxGetDimensions(prhs[3]);
	if (nx != dim_array[0] || ny != dim_array[1]) {
		mexErrMsgTxt("The size of the phase images are not equal.");
	}
	dim_array = mxGetDimensions(prhs[4]);
	if (nx != dim_array[0] || ny != dim_array[1]) {
		mexErrMsgTxt("The size of the phase images are not equal.");
	}
	dim_array = mxGetDimensions(prhs[5]);
	if (nx != dim_array[0] || ny != dim_array[1]) {
		mexErrMsgTxt("The size of the weight1 are not right. It must be equal to the size of phase image.");
	}
	dim_array = mxGetDimensions(prhs[6]);
	if (nx != dim_array[0] || ny != dim_array[1]) {
		mexErrMsgTxt("The size of the weight2 are not right. It must be equal to the size of phase image");
	}

	// allocate and copy the input HARP maps
	phase1_t0 = pgDmatrix(0, ny-1, 0, nx-1);
	phase1_t1 = pgDmatrix(0, ny-1, 0, nx-1);
	phase2_t0 = pgDmatrix(0, ny-1, 0, nx-1);
	phase2_t1 = pgDmatrix(0, ny-1, 0, nx-1);
	weight1	  =	pgDmatrix(0, ny-1, 0, nx-1);
	weight2	  =	pgDmatrix(0, ny-1, 0, nx-1);

	// copy the input HARP images to the memory
	pd = (double*)mxGetPr(prhs[1]);
	memcpy(&(phase1_t0[0][0]), pd, nx*ny*sizeof(double));

	pd = (double*)mxGetPr(prhs[2]);
	memcpy(&(phase1_t1[0][0]), pd, nx*ny*sizeof(double));

	pd = (double*)mxGetPr(prhs[3]);
	memcpy(&(phase2_t0[0][0]), pd, nx*ny*sizeof(double));

	pd = (double*)mxGetPr(prhs[4]);
	memcpy(&(phase2_t1[0][0]), pd, nx*ny*sizeof(double));

	pd = (double*)mxGetPr(prhs[5]);
	memcpy(&(weight1[0][0]), pd, nx*ny*sizeof(double));
	
	pd = (double*)mxGetPr(prhs[6]);
	memcpy(&(weight2[0][0]), pd, nx*ny*sizeof(double));
	
	//normalize the weight
	maxWeight1 = 0; maxWeight2 = 0;
	for (j=0; j<ny; j++)
	{
		for (i=0; i<nx; i++)
		{
			if (fabs(weight1[j][i]) > maxWeight1)
				maxWeight1 = fabs(weight1[j][i]);
			if (fabs(weight2[j][i]) > maxWeight2)
				maxWeight2 = fabs(weight2[j][i]);
		}
	}
	for (j=0; j<ny; j++)
	{
		for (i=0; i<nx; i++)
		{
			weight1[j][i] = fabs(weight1[j][i])/maxWeight1;
			weight2[j][i] = fabs(weight2[j][i])/maxWeight2;
		}
	}
	    
	// the sixth input: the position of the seed at the second time frame
	if (nrhs == 7)
	{
		// track the seed
		ptPhase[0] = BilinearInterpolationWrap(phase1_t0, nx, ny, startPoint[0], startPoint[1]);
		ptPhase[1] = BilinearInterpolationWrap(phase2_t0, nx, ny, startPoint[0], startPoint[1]);
		ret = TrackPoint(phase1_t1, phase2_t1, nx, ny, startPoint, ptPhase, trackedPt, motion);
		
	} else if (nrhs == 8) {
	
		pd = (double*)mxGetPr(prhs[7]);
		dim_array = mxGetDimensions(prhs[7]);
		if ( (dim_array[0] == 1 && dim_array[1] == 2) || 
			 (dim_array[0] == 2 && dim_array[1] == 1))	// means if the 8th input is the tracked seed
		{
			trackedPt[0] = pd[0]-1;
			trackedPt[1] = pd[1]-1;
			
		} else {			// means if the last input is not 2 by 1 vector (tracked seed)
				// track the seed
			ptPhase[0] = BilinearInterpolationWrap(phase1_t0, nx, ny, startPoint[0], startPoint[1]);
			ptPhase[1] = BilinearInterpolationWrap(phase2_t0, nx, ny, startPoint[0], startPoint[1]);
			ret = TrackPoint(phase1_t1, phase2_t1, nx, ny, startPoint, ptPhase, trackedPt, motion);
		}
	}
    
    // need allocate memory for those input variables. 
    motion_x = pgDmatrix(0, ny-1, 0, nx-1);
    motion_y = pgDmatrix(0, ny-1, 0, nx-1);
    trackIndex = pgDmatrix(0, ny-1, 0, nx-1);
	savedWeight = pgDmatrix(0, ny-1, 0, nx-1);
	savedWeight1 = pgDmatrix(0, ny-1, 0, nx-1);    
	savedWeight2 = pgDmatrix(0, ny-1, 0, nx-1);
	
	savedWeight[0][0] = dim_array[0];
	savedWeight[0][1] = dim_array[1];

    TrackingRefinement_FastMarching(phase1_t0, phase1_t1, phase2_t0, phase2_t1, weight1, weight2, 
									nx, ny, startPoint, trackedPt, motion_x, motion_y, trackIndex, savedWeight,savedWeight1,savedWeight2);
    //rgIndex = IdentifyGoodRegion(motion_x, motion_y, trackIndex, startPoint, nx, ny);
		
	plhs[0] = mxCreateDoubleMatrix(nx, ny, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nx, ny, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nx, ny, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(nx, ny, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(nx, ny, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(nx, ny, mxREAL);
    
    memcpy((char*)mxGetData(plhs[0]), &(motion_x[0][0]), nx*ny*sizeof(double));
    memcpy((char*)mxGetData(plhs[1]), &(motion_y[0][0]), nx*ny*sizeof(double));
    memcpy((char*)mxGetData(plhs[2]), &(trackIndex[0][0]), nx*ny*sizeof(double));

	memcpy((char*)mxGetData(plhs[3]), &(savedWeight[0][0]), nx*ny*sizeof(double));
	memcpy((char*)mxGetData(plhs[4]), &(savedWeight1[0][0]), nx*ny*sizeof(double));
	memcpy((char*)mxGetData(plhs[5]), &(savedWeight2[0][0]), nx*ny*sizeof(double));
    //memcpy((char*)mxGetData(plhs[2]), &(rgIndex[0][0]), nx*ny*sizeof(double));
	
	
    pgFreeDmatrix(savedWeight,0, ny-1, 0, nx-1);
	pgFreeDmatrix(savedWeight1,0, ny-1, 0, nx-1);
	pgFreeDmatrix(savedWeight2,0, ny-1, 0, nx-1);
    
    pgFreeDmatrix(motion_x,0, ny-1, 0, nx-1);
    pgFreeDmatrix(motion_y,0, ny-1, 0, nx-1);
    pgFreeDmatrix(trackIndex,0, ny-1, 0, nx-1);
    
    //pgFreeDmatrix(rgIndex, 0, ny-1, 0, nx-1);
    

	// delete phase1_t0, phase1_t1, phase2_t0, phase2_t1
	pgFreeDmatrix(phase1_t0,0, ny-1, 0, nx-1);
	pgFreeDmatrix(phase1_t1,0, ny-1, 0, nx-1);
	pgFreeDmatrix(phase2_t0,0, ny-1, 0, nx-1);
	pgFreeDmatrix(phase2_t1,0, ny-1, 0, nx-1);
	pgFreeDmatrix(weight1,0, ny-1, 0, nx-1);
	pgFreeDmatrix(weight2,0, ny-1, 0, nx-1);
	
	
    
	return;


}

double BilinearInterpolation(double **mat, int XO, int YO, double x, double y)
{
    int i0, j0, i1, j1;
    double dx, dy, hx, hy;
    double result;

    if(x < 0 || x > (XO-1) || y < 0 || y > (YO-1) )
    {
        result = 0;
    }
    else
    {
        j1 = (int)ceil(x);
        i1 = (int)ceil(y);
        j0 = (int)floor(x);
        i0 = (int)floor(y);
        dx = x - (float)j0;
        dy = y - (float)i0;
        /* Introduce more variables to reduce computation */
        hx = 1.0 - dx;
        hy = 1.0 - dy;

        /* Trilinear Interpolate of vx */
        /* 
            Vxy = V00(1-x)(1-y) + 
            V10x(1-y) +
            V01(1-x)y +
            V11xy;
        */
        result = mat[i0][j0]*(hx)*(hy) +
            mat[i0][j1]*dx*(hy) +
            mat[i1][j0]*(hx)*dy +
            mat[i1][j1]*dx*dy;
       //printf("oldV[%d][%d] = %d result = %f, x=%f, y=%f\n", i0, j0, oldV[i0][j0], result, x, y);
    }
      
    return result;
}


// *********** Bilinear interpolation of the phase

double BilinearInterpolationWrap(double **mat, int sx, int sy, double x, double y)
{
	int i0, j0, i1, j1;
	double dx, dy, hx,hy;
    
    double ta, tb, tc, val;

	j1 = (int)ceil(x);
	i1 = (int)ceil(y);
	j0 = (int)floor(x);
	i0 = (int)floor(y);
    
	if (j0<0 || j1>(sx-1) || i0<0 || i1>(sy-1))
	{
		return 0;
	} else {
		dx = x-j0;
		dy = y-i0;
		hx = 1.0-dx;
		hy = 1.0-dy;
		//return (mat[i0][j0]*hx + mat[i0][j1]*dx)*hy + (mat[i1][j0]*hx + mat[i1][j1]*dx)*dy;
		
		tb = WRAP(mat[i1][j0] - mat[i0][j0]);
		tc = WRAP(mat[i0][j1] - mat[i0][j0]);
		ta = -tb - tc + WRAP(mat[i1][j1] - mat[i0][j0]);
		val = WRAP(ta*dx*dy + tb*dy + tc*dx + mat[i0][j0]);

		return val;

	}

}

double BilinearInterpolationWrap2(double **mat, int sx, int sy, double x, double y, double *grad)
{
	int i0, j0, i1, j1;
	double dx, dy, hx,hy;
    double ta, tb, tc, val;

	j0 = (int)floor(x);
	if (j0 == sx-1)
		j0 = j0-1;
	
	i0 = (int)floor(y);
	if (i0 == sy-1)
		i0 = i0-1;
		
	j1 = j0+1;
	i1 = i0+1;

	if (j0<0 || j1>(sx-1) || i0<0 || i1>(sy-1))
	{
		return 0;
	} else {
		dx = x-j0;
		dy = y-i0;
		hx = 1.0-dx;
		hy = 1.0-dy;
		//return (mat[i0][j0]*hx + mat[i0][j1]*dx)*hy + (mat[i1][j0]*hx + mat[i1][j1]*dx)*dy;
		
		tb = WRAP(mat[i1][j0] - mat[i0][j0]);
		tc = WRAP(mat[i0][j1] - mat[i0][j0]);
		ta = -tb - tc + WRAP(mat[i1][j1] - mat[i0][j0]);
		val = WRAP(ta*dx*dy + tb*dy + tc*dx + mat[i0][j0]);
		
        grad[0] = tc + ta*dy;
		grad[1] = tb + ta*dx;

		return val;

	}

}


int TrackPoint(double **phase1_t1, double **phase2_t1, 
				  int nx, int ny,
				  double* point2D, 
				  double* currentPhase,
				  double* trackedPt, 
                  double* motion
				  )
{

	int maxIter = 15;

	int n;
    
    double phase[2];
    double phaseDiff[2];

	double grad0[2], grad1[2];
    
    double det;

	if (point2D[0]<0 || point2D[0]>nx-1 || point2D[1]<0 || point2D[1]>ny-1)
		return -1;
	

	//double trackedPt[2];
	trackedPt[0] = point2D[0];
	trackedPt[1] = point2D[1];

	for (n = 0; n<maxIter; n++) 
	{
		// bilinear interpolation with wrapping
		phase[0] = BilinearInterpolationWrap2(phase1_t1, nx, ny, trackedPt[0], trackedPt[1], grad0);
		phase[1] = BilinearInterpolationWrap2(phase2_t1, nx, ny, trackedPt[0], trackedPt[1], grad1);

		phaseDiff[0] = WRAP(phase[0] - currentPhase[0]);
		phaseDiff[1] = WRAP(phase[1] - currentPhase[1]);
		
		det = grad0[0]*grad1[1] - grad0[1]*grad1[0];
		motion[0] = -( grad1[1]*phaseDiff[0] - grad0[1]*phaseDiff[1]) / det;
		motion[1] = -(-grad1[0]*phaseDiff[0] + grad0[0]*phaseDiff[1]) / det;
		if (fabs(motion[0]) > 1)
			motion[0] = SIGN(motion[0]);
		if (fabs(motion[1]) > 1)
			motion[1] = SIGN(motion[1]);

		trackedPt[0] = trackedPt[0] + motion[0];
		trackedPt[1] = trackedPt[1] + motion[1];

		if (fabs(phaseDiff[0]) + fabs(phaseDiff[1]) <=0.0001)
			break;
	}
    motion[0] = trackedPt[0] - point2D[0];
	motion[1] = trackedPt[1] - point2D[1];
	if (n == maxIter)
	{
		//trackedPt[0] = point2D[0];
		//trackedPt[1] = point2D[1];
        //motion[0] = -1000;
        //motion[1] = -1000;
		return -1;
	} else {
        
		return 0;
	}


}



void TrackingRefinement_FastMarching(double **phase1_t0, double **phase1_t1, double **phase2_t0, double **phase2_t1,
						double **weight1, double **weight2,
                        int nx, int ny,
                        double *startPt,
						double *trackedPt,
                        double **motion_x,
                        double **motion_y, 
                        double **trackIndex,
						double **savedWeight, double **savedWeight1,double **savedWeight2
                        )
                        
{
    double currentPhase[2];
    int istartPt[2];
    double currentPt[2];
	double estPt[2];
	double closestPt[2];
    
    double ptPhase[2];
    //double trackedPt[2];
    double ptMotion[2];
    SSDL_node *startNode, *currentNode, *nextNode, *tempNode, *newNode, *node;
	
	SSDL_node ***nodeMatrix;  // each position connect to a broundary point nood

    int neighbor_x[8];
    int neighbor_y[8];
    double neighbor_val[8];
    int index[8];    // used for sorting
    
    int no_neighbor;
    
    //double ** delete
    
    int val, ret;
    
    int i, j, x, y, n;
	int ix, iy;
    double oldCost, newCost;
    int num_Neighbor;
    
    double temp_x, temp_y, temp;
    double temp_phase1, temp_phase2;
	int neighbors[8][2];
    
    //int **trackIndex;
    int linkLen;
	int count;
	char filename[128];
	char tempStr[10];
	
	double closestMotion[2];
	
	int l0, l1, l2, l3;
	
	double tempWeight, tempWeight1, tempWeight2;
	double addedCost;
	
	nodeMatrix = (SSDL_node ***) pgULmatrix(0,ny-1,0,nx-1);

	// the relative position of four neighbors
	neighbors[0][0] = -1; neighbors[0][1] =  0;
	neighbors[1][0] =  1; neighbors[1][1] =  0;
	neighbors[2][0] =  0; neighbors[2][1] =  1;
	neighbors[3][0] =  0; neighbors[3][1] = -1;
    neighbors[4][0] = -1; neighbors[4][1] = -1;
	neighbors[5][0] = -1; neighbors[5][1] =  1;
	neighbors[6][0] =  1; neighbors[6][1] = -1;
	neighbors[7][0] =  1; neighbors[7][1] =  1;
    
	// initialize matrix trackIndex    
	// index matrix: 
    //          1 means tracked points, 
    //          2 means points in the SSL, 
    //          0 means not-tracked points, and
    //          -1 means points cannot be tracked correct (due to phase inconsistency)
    for (j=0; j<ny; j++)
    {
        for (i=0; i<nx; i++)
        {
            trackIndex[j][i] = 0;
			motion_x[j][i] = 0;
			motion_y[j][i] = 0;
			nodeMatrix[j][i] = NULL;
        }
    }
    
	// **********   Step 1: Initialization ****************//
	
	// 1.1 track the rounded seed point

    currentPt[0] = (int)floor(startPt[0]);
    currentPt[1] = (int)floor(startPt[1]);
    
    ptPhase[0] = BilinearInterpolationWrap(phase1_t0, nx, ny, currentPt[0], currentPt[1]);
	ptPhase[1] = BilinearInterpolationWrap(phase2_t0, nx, ny, currentPt[0], currentPt[1]);
	
	estPt[0] = currentPt[0] + trackedPt[0] - startPt[0];
	estPt[1] = currentPt[1] + trackedPt[1] - startPt[1];
	
	// assume the point is already tracked and stored in trackedPt, so needn't track the seed again.
	
	ptMotion[0] = trackedPt[0] - startPt[0];
	ptMotion[1] = trackedPt[1] - startPt[1];
	// here assume the point can always be tracked correctly. Need take more care later 
	ret = TrackPoint(phase1_t1, phase2_t1, nx, ny, estPt, ptPhase, trackedPt, ptMotion);
	if (val == -1) // point track failed
        return;
	
	motion_x[(int)currentPt[1]][(int)currentPt[0]] = trackedPt[0] - currentPt[0];
	motion_y[(int)currentPt[1]][(int)currentPt[0]] = trackedPt[1] - currentPt[1];
	trackIndex[(int)currentPt[1]][(int)currentPt[0]] = 1;


	// 1.2, put its neighbors into the edge list, and the SSDL
    num_Neighbor = 4;       // 4-neighbor
	//num_Neighbor = 8;       // 8-neighbor
	
	
	// put the seed into the SSL
	startNode = NULL;
	linkLen = 0;
	
	nextNode = (SSDL_node *)malloc(sizeof(SSDL_node));
	nextNode->x					= (int)(currentPt[0]);
	nextNode->y					= (int)(currentPt[1]);
	nextNode->accumulatedCost	= 0;
	nextNode->closest_x			= (int)(currentPt[0]);
	nextNode->closest_y			= (int)(currentPt[1]);
	nextNode->prev				= NULL;
	nextNode->next				= NULL;

    startNode = InsertSSDL(startNode, nextNode);
	linkLen ++;
	nodeMatrix[nextNode->y][nextNode->x] = nextNode;
	trackIndex[nextNode->y][nextNode->x] = 2;
	
	// index matrix trackIndex: 
	//          1 means tracked points, 
	//          2 means boundary points, 
	//          0 means not-visited points, and
	//         -1 means points cannot be tracked correctly (due to phase inconsistency)

    
    // **************** step 2: main body of region growing ************//
    
	

	count = 0; // used for saving intermediate result into file
    while (startNode != NULL) // or while (linkLen > 0)
//for (n=0; n<60; n++)
    {
		// 2.2: track the point at the top of the list
		currentNode = startNode;
        startNode = startNode->next;
		if (startNode != NULL)
			startNode->prev = NULL;
        
		closestPt[0] = currentNode->closest_x;
		closestPt[1] = currentNode->closest_y;
		
        currentPt[0] = currentNode->x + motion_x[currentNode->closest_y][currentNode->closest_x]; // start tracking position
        currentPt[1] = currentNode->y + motion_y[currentNode->closest_y][currentNode->closest_x];
        ptPhase[0] = phase1_t0[currentNode->y][currentNode->x];  // phase value in initial time frame
        ptPhase[1] = phase2_t0[currentNode->y][currentNode->x];

		// track the current point, starting by assuming same amount of motion as the closest point
        ret = TrackPoint(phase1_t1, phase2_t1, nx, ny, currentPt, ptPhase, trackedPt, ptMotion);
		
//		motion_x[0][0] = currentNode->accumulatedCost;

/*		if (currentNode->x == 45 && currentNode->y ==79)
		{
			savedWeight[0][0] = currentNode->x+1;
			savedWeight[0][1] = currentNode->y+1;
			savedWeight[0][2] = currentNode->closest_x+1;
			savedWeight[0][3] = currentNode->closest_y+1;
			savedWeight[0][4] = currentPt[0]+1;
			savedWeight[0][5] = currentPt[1]+1;
			savedWeight[0][6] = ptPhase[0];
			savedWeight[0][7] = ptPhase[1];
			savedWeight[0][8] = trackedPt[0]+1;
			savedWeight[0][9] = trackedPt[1]+1;

		}//*/
			
        if (ret == -1 || (sqrt(ptMotion[0]*ptMotion[0] + ptMotion[1]*ptMotion[1])) > 1)
		//if (ret == -1)
        {
            trackIndex[currentNode->y][currentNode->x] = -1; // means the point cannot be correctly tracked
			nodeMatrix[currentNode->y][currentNode->x] = NULL;
			free(currentNode);
			linkLen --;

            continue;
        }
		// store tracked motion
		motion_x[currentNode->y][currentNode->x] = trackedPt[0] - currentNode->x;
        motion_y[currentNode->y][currentNode->x] = trackedPt[1] - currentNode->y;;
        // mark this point as tracked.
        trackIndex[currentNode->y][currentNode->x] = 1;
		count ++;
/*		
	if (floor(count/50)*50 == count)
		{
			strcpy(filename, "temp_FM");
			sprintf(tempStr, "%03d", (int)floor(count/50));
			strcat(filename, tempStr);
			strcat(filename, ".dat");
			pgWriteDouble(&(trackIndex[0][0]), ny*nx, filename);

		}//*/
		
		
		// 2.2: check the neighbors, and update their cost if necessary
		
		for (i=0; i<num_Neighbor; i++)
		{
			ix = currentNode->x + neighbors[i][0];
			iy = currentNode->y + neighbors[i][1];
			
			if (  ix < 0 || ix > nx-1 || iy < 0 || iy > ny-1)
				continue;
			
				
			if (trackIndex[iy][ix] == 1 || trackIndex[iy][ix] == -1) // means already tracked
				continue;

			temp_x = (double)ix + motion_x[currentNode->y][currentNode->x];		// assume same motion as the currentNode
			temp_y = (double)iy + motion_y[currentNode->y][currentNode->x];
			if (temp_x<0 || temp_x>nx-1 || temp_y<0 || temp_y>ny-1)
				continue;
			
			temp_phase1 = BilinearInterpolationWrap(phase1_t1, nx, ny, temp_x, temp_y);
			temp_phase2 = BilinearInterpolationWrap(phase2_t1, nx, ny, temp_x, temp_y);
				
			tempWeight1 = weight1[iy][ix];
			tempWeight2 = BilinearInterpolation(weight2, nx, ny, temp_x, temp_y);
			//tempWeight = 1/(0.0001+(tempWeight1 + tempWeight2));
			tempWeight = 1/(0.0001 + tempWeight1*tempWeight2);
			
			//for debug
			if (iy != 0) {
				savedWeight[iy][ix] = tempWeight;
				savedWeight1[iy][ix] = tempWeight1;
				savedWeight2[iy][ix] = tempWeight2;
			}

			addedCost = fabs(WRAP(temp_phase1 - phase1_t0[iy][ix])) + 
						fabs(WRAP(temp_phase2 - phase2_t0[iy][ix]));
						//pow(fabs(WRAP(temp_phase1 - phase1_t0[iy][ix])),2) + 
						//pow(fabs(WRAP(temp_phase2 - phase2_t0[iy][ix])),2);
			addedCost *= tempWeight;
			
			if (i<4) // means the 4 closer neighbor
			{
				newCost = currentNode->accumulatedCost + addedCost;				
			} else {
				newCost = currentNode->accumulatedCost + sqrt(2.0)*addedCost;
			}

			if ( trackIndex[iy][ix] == 2) // if the neighbor point is in the broundary list,update its cost
			{
				node = nodeMatrix[iy][ix];
				// if newcost is more, do nothing
				if (newCost >= node->accumulatedCost)
					continue;
			
				// otherwise, update its accumulatedCost
				
				//l1 = SSDLLength(startNode);
				// if the new cost is smaller, update the node
				if (node->prev != NULL) 
					(node->prev)->next = node->next;
				
				if (node->next != NULL)
					(node->next)->prev = node->prev;
					
				if (node == startNode)
					startNode = node->next;
					
/*					l2 = SSDLLength(startNode);
				if (l2 != l1-1)
				{
					trackIndex[0][0] = linkLen;
					trackIndex[0][1] = l1;
					trackIndex[0][2] = l2;
					if (startNode == NULL)
						trackIndex[0][3] = -100;
					if (node->prev == NULL)
						trackIndex[0][3] = -200;
					if (node->next == NULL)
						trackIndex[0][3] = -300;
					return;
				}
*/
				node->accumulatedCost = newCost;
				
				node->closest_x = currentNode->x;
				node->closest_y = currentNode->y;
				node->prev = NULL;
				node->next = NULL;
				
				l1 = SSDLLength(startNode);
				startNode = InsertSSDL(startNode, node);
						

			} else if (trackIndex[iy][ix] == 0) { // if point never accessed, then add it to the list 
				
				newNode = (SSDL_node *)malloc(sizeof(SSDL_node));
				newNode->x = ix;
				newNode->y = iy;
				
				newNode->accumulatedCost = newCost;
			
				newNode->closest_x = currentNode->x;
				newNode->closest_y = currentNode->y;
				newNode->prev = NULL;
				newNode->next = NULL;
				
				startNode = InsertSSDL(startNode, newNode);
				
				trackIndex[newNode->y][newNode->x] = 2;
				
				nodeMatrix[newNode->y][newNode->x] = newNode;
				linkLen ++;
				
			}
		} // end for neighbors
		
		
		nodeMatrix[currentNode->y][currentNode->x] = NULL;
        free(currentNode);
		linkLen --;
        
    }// end while loop
/*
	strcpy(filename, "temp_FM");
	sprintf(tempStr, "%03d", (int)floor(count/50)+1);
	strcat(filename, tempStr);
	strcat(filename, ".dat");
	pgWriteDouble(&(trackIndex[0][0]), ny*nx, filename);
    //*/
	
	trackIndex[0][0] = linkLen;
    
    // free memory
    currentNode = startNode;
    while ( currentNode != NULL)
    {
        nextNode = currentNode->next;
        free(currentNode);
        currentNode = nextNode;
    }
	//*/
	pgFreeULmatrix((unsigned long long **)nodeMatrix, 0,ny-1,0,nx-1);
    return;   
}
 
// find the link length

int SSDLLength(SSDL_node *start)
{
	int linkLen;

	linkLen = 0;
	while (start != NULL)
	{
		linkLen ++;
		start = start->next;
	}
	
	return linkLen;
}

// ** insert a new node into the sequentially sorted link
SSDL_node* InsertSSDL(SSDL_node *start, SSDL_node *newNode)
{
    int n;
    SSDL_node *next, *prev, *current;
	SSDL_node *newStart;
	
	if (start == NULL)
	{
		newNode->next = NULL;
		newNode->prev = NULL;
		return newNode;
	}
    
	newStart = start;
    if (start->accumulatedCost > newNode->accumulatedCost)
    {
		newNode->prev = NULL;
        newNode->next = start;
		start->prev = newNode;
        return newNode;
    } else {
        current = start->next;
        prev = start;
        while (current != NULL)
        {
            if (current->accumulatedCost > newNode->accumulatedCost)
                break;
            prev = current;
            current = prev->next;
        }
        // insert newNode
        prev->next = newNode;
		newNode->prev = prev;
        newNode->next = current;
		if (current != NULL) 
		{
			current->prev = newNode;
		}
    }
    return newStart;
}

// quick sort
void QKSort3(int n, double* arr, int *brr, int *crr)
{
    int l=0, jstack=0,ir,iq,i,j,k;
    int istack[50];
    long int fx=0L;
    double a, temp;
    int b, c;

    ir=n-1;
    for (;;) {

        if (ir-l < 7) {
            for (j=l+1;j<=ir;j++) {
	            a = arr[j];
	            b = brr[j];
                c = crr[j];
	            for (i=j-1;arr[i]>a && i>=0;i--){ 
        	        arr[i+1] = arr[i];
	                brr[i+1] = brr[i];
                    crr[i+1] = crr[i];
	            }
	            arr[i+1] = a;
	            brr[i+1] = b;
                crr[i+1] = c;
            }
            if (jstack == 0) return;
            ir=istack[jstack--];
            l=istack[jstack--];
        } 
         else {
              k = (l+ir) >> 1;
              swap(double, arr[k], arr[ir]);
              swap(int, brr[k], brr[ir]);
              swap(int, crr[k], crr[ir]);
              if (arr[l] > arr[ir]){
                  swap(double, arr[l], arr[ir]);
                  swap(int, brr[l], brr[ir]);
                  swap(int, crr[l], crr[ir]);
              }
              if (arr[l+1] > arr[ir]){
                  swap(double, arr[l+1], arr[ir]);
                  swap(int, brr[l+1], brr[ir]);
                  swap(int, crr[l+1], crr[ir]);
              }
              if (arr[l] > arr[l+1]){
                  swap(double, arr[l+1], arr[l]);
                  swap(int, brr[l+1], brr[l]);
                  swap(int, crr[l+1], crr[l]);
              }
              i = l+1;
              j = ir;
              a = arr[l+1];
              b = brr[l+1];
              c = crr[l+1];
              for (;;) {
                  while (j >= 0 && a < arr[j]) 
                      j--;
                  if (j <= i) {
                      arr[i]=a;
                      brr[i]=b;
                      crr[i]=c;
                      break;
                  }
                  arr[i]=arr[j];
                  brr[i]=brr[j];
                  crr[i]=crr[j];
                  i++;
                  while (a > arr[i] && i < n) 
                      i++;
                  if (j <= i) {
                      arr[(i=j)]=a;
                      brr[i] = b;
                      crr[i] = c;
                      break;
                  }
                  arr[j]=arr[i];
                  brr[j]=brr[i];	
                  crr[j]=crr[i];
                  j--;
              }
              if (ir-i >= i-l) {
                  istack[++jstack]=i+1;
                  istack[++jstack]=ir;
                  ir=i-1;
              } else {
                  istack[++jstack]=l;
                  istack[++jstack]=i-1;
                  l=i+1;
              }
         }
    }
}
    
    

void QKSort2(int n, double* arr, int *brr)
{
    int l=0, jstack=0,ir,iq,i,j,k;
    int istack[50];
    long int fx=0L;
    double a, temp;
    int b;

    ir=n-1;
    for (;;) {

        if (ir-l < 7) {
            for (j=l+1;j<=ir;j++) {
	            a = arr[j];
	            b = brr[j];
	            for (i=j-1;arr[i]>a && i>=0;i--){ 
        	        arr[i+1] = arr[i];
	                brr[i+1] = brr[i];
	            }
	            arr[i+1] = a;
	            brr[i+1] = b;
            }
            if (jstack == 0) return;
            ir=istack[jstack--];
            l=istack[jstack--];
        } 
         else {
              k = (l+ir) >> 1;
              swap(double, arr[k], arr[ir]);
              swap(int, brr[k], brr[ir]);
              if (arr[l] > arr[ir]){
                  swap(double, arr[l], arr[ir]);
                  swap(int, brr[l], brr[ir]);
              }
              if (arr[l+1] > arr[ir]){
                  swap(double, arr[l+1], arr[ir]);
                  swap(int, brr[l+1], brr[ir]);
              }
              if (arr[l] > arr[l+1]){
                  swap(double, arr[l+1], arr[l]);
                  swap(int, brr[l+1], brr[l]);
              }
              i = l+1;
              j = ir;
              a = arr[l+1];
              b = brr[l+1];
              for (;;) {
                  while (j >= 0 && a < arr[j]) 
                      j--;
                  if (j <= i) {
                      arr[i]=a;
                      brr[i]=b;
                      break;
                  }
                  arr[i]=arr[j];
                  brr[i]=brr[j];
                  i++;
                  while (a > arr[i] && i < n) 
                      i++;
                  if (j <= i) {
                      arr[(i=j)]=a;
                      brr[i] = b;
                      break;
                  }
                  arr[j]=arr[i];
                  brr[j]=brr[i];	
                  j--;
              }
              if (ir-i >= i-l) {
                  istack[++jstack]=i+1;
                  istack[++jstack]=ir;
                  ir=i-1;
              } else {
                  istack[++jstack]=l;
                  istack[++jstack]=i-1;
                  l=i+1;
              }
         }
    }
}


double **IdentifyGoodRegion(double **motion_x, double **motion_y, double **trackIndex, double *startPt, int nx, int ny)
{
    int i, j, x, y;
    RG_node *start, *prev, *next, *newNode, *current;
    
    double** rgIndex;
    int neighbor[4][2];
    int linkLen;
    int startx, starty;
    
    // the threshold for region growing. When motion between neighboring points are larger than thre,
    // they are thought to be in different regions.
   
    double thre;
    thre = 0.8;
    
    startx = (int)floor(startPt[0]);
    starty = (int)floor(startPt[1]);
    
    neighbor[0][0] = -1; neighbor[0][1] = 0;
    neighbor[1][0] =  1; neighbor[1][1] = 0;
    neighbor[2][0] =  0; neighbor[2][1] = -1;
    neighbor[3][0] =  0; neighbor[3][1] = 1;
    
    rgIndex = pgDmatrix(0, ny-1, 0, nx-1);
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            rgIndex[j][i] = 0;
        }
    }
    
    start = (RG_node*)malloc(sizeof(RG_node));
    start->x = startx;
    start->y = starty;
    start->next = NULL;
	rgIndex[start->y][start->x] = 1;
    
    linkLen = 1;
    
    while(start != NULL)
    {
        current = start;
        start = start->next;
        linkLen --;
        
        for (i = 0; i<4; i++) 
        {
            x = current->x + neighbor[i][0];
            y = current->y + neighbor[i][1];
            
			if (   x<0 || x>nx-1 || y<0 || y>ny-1
				|| rgIndex[y][x] != 0
				|| trackIndex[y][x] != 1
			   )
			{

			} else {
				if (fabs(motion_x[current->y][current->x]-motion_x[y][x]) < thre &&
					fabs(motion_y[current->y][current->x]-motion_y[y][x]) < thre)
				{
					newNode = (RG_node*)malloc(sizeof(RG_node));
	                newNode->x = x;
		            newNode->y = y;
			        newNode->next = start;
				    start = newNode;
                    linkLen ++;

					rgIndex[y][x] = 1;
				} else {
					rgIndex[y][x] = -1;
				}
		    }
		}
        
        free(current);
    }
    
    return rgIndex;
}



int pgWriteDouble(double *ptr, int size, const char *filename)
{
  FILE *fp;		
  int numwritten;

  if (filename == NULL) 
    fp = stdout;
  else {
    if ((fp=fopen(filename,"wb")) == NULL) 
      return -1;
  }
  numwritten = (int) fwrite((char *)ptr,sizeof(double),size,fp);
  fclose(fp);

  if (numwritten != size) 
    return -1;
  else
    return 0;
}

/* Allocates a matrix of doubles with range [nrl..nrh][ncl..nch] */
double **pgDmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  
  unsigned long long bufsize,bufptr;
  double **m;

  bufsize = (nrh-nrl+1)*sizeof(double*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(double);

  m=(double **) malloc(bufsize);
  if (!m) {
	  return NULL;
  }
  m -= nrl;

  bufptr = ((unsigned long long) (m+nrl)) + (nrh-nrl+1)*sizeof(double*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((double *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(double)));
    m[j] -= ncl;
  }

  return m;
}



/* Frees a matrix allocated by dmatrix */
void pgFreeDmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}

/* Allocates a matrix of doubles with range [nrl..nrh][ncl..nch] */
int **pgImatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  
  unsigned long long bufsize,bufptr;
  int **m;

  bufsize = (nrh-nrl+1)*sizeof(int*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(int);

  m=(int **) malloc(bufsize);
  if (!m) {
	  return NULL;
  }
  m -= nrl;

  bufptr = ((unsigned long long) (m+nrl)) + (nrh-nrl+1)*sizeof(int*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((int *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(int)));
    m[j] -= ncl;
  }

  return m;
}



/* Frees a matrix allocated by dmatrix */
void pgFreeImatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}

/* Allocates a matrix of unsigned chars with range [nrl..nrh][ncl..nch] */
char **pgBmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  unsigned long long bufsize,bufptr;
  char **m;

  bufsize = (nrh-nrl+1)*sizeof(char*)
            + (nrh-nrl+1)*(nch-ncl+1)*sizeof(char);

  m=(char **) malloc(bufsize);
  if (!m) return NULL;
  m -= nrl;

  bufptr = ((unsigned long long) (m+nrl)) + (nrh-nrl+1)*sizeof(char*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((char *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(char)));
    m[j] -= ncl;
  }

  return m;
}




/* Frees a matrix allocated by bmatrix */
void pgFreeBmatrix(char **m,int nrl,int nrh,int ncl,int nch)
{
        free((char*) (m+nrl));
}

/* Allocates a matrix of doubles with range [nrl..nrh][ncl..nch] */
unsigned long long **pgULmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  
  unsigned long long bufsize,bufptr;
  unsigned long long **m;

  bufsize = (nrh-nrl+1)*sizeof(unsigned long long*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(unsigned long long);

  m=(unsigned long long **) malloc(bufsize);
  if (!m) {
	  return NULL;
  }
  m -= nrl;

  bufptr = ((unsigned long long) (m+nrl)) + (nrh-nrl+1)*sizeof(unsigned long long*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((unsigned long long *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(unsigned long long)));
    m[j] -= ncl;
  }

  return m;
}



/* Frees a matrix allocated by dmatrix */
void pgFreeULmatrix(unsigned long long **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}




