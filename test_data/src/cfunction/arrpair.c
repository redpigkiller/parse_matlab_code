#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/*
    This function performs pairing between two arrays with same size, 
    it minimizes the "circular" distance and returns the paired second array.
    Make sure the input range is [-cirbnd, cirbnd], otherwise,
    the distance between two points may be negative. (cirbnd = inf by default)

    Algorithm: Hungarian Algorithm / Munkres algorithm

	Usage:
		[arr_hat, indices, cost] = arrpair(arr_ref, arr_hat, cir_bnd, method);

		+ input args:
			arr_ref: 	input vector of existing angle
			arr_hat:	input vector of measured angle
			cir_bnd:	the angle is periodic. By default, cir_bnd = inf, i.e., the value of angle is in [-0.5, 0.5). For cir_bnd < 0, then there is no circular bound.
			method:		0 for assignmentoptimal (default)
						1 for assignmentsuboptimal1
						2 for assignmentsuboptimal2
		+ output args:
			arr_hat:		output vector of measured arr_hat paired with arr_ref such that 1-norm of arr_hat-arr_ref is minimized
			indices:		indices mapping from arr_hat to arr_hat, i.e., arr_hat = arr_hat(indices) in MATLAB
			cost:			1-norm of arr_hat-arr_ref
    
    ref:
        https://github.com/mcximing/hungarian-algorithm-cpp
        https://github.com/phoemur/hungarian_algorithm
        https://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem
*/
/*
ASSIGNMENTOPTIMAL    Compute optimal assignment by Munkres algorithm
   ASSIGNMENTOPTIMAL(DISTMATRIX) computes the optimal assignment (minimum
   overall costs) for the given rectangular distance or cost matrix, for
   example the assignment of tracks (in rows) to observations (in
   columns). The result is a column vector containing the assigned column
   number in each row (or 0 if no assignment could be done).

   [ASSIGNMENT, COST] = ASSIGNMENTOPTIMAL(DISTMATRIX) returns the
   assignment vector and the overall cost.

   The distance matrix may contain infinite values (forbidden
   assignments). Internally, the infinite values are set to a very large
   finite number, so that the Munkres algorithm itself works on
   finite-number matrices. Before returning the assignment, all
   assignments with infinite distance are deleted (i.e. set to zero).

   A description of Munkres algorithm (also called Hungarian algorithm)
   can easily be found on the web.

   <a href="assignment.html">assignment.html</a>  <a href="http://www.mathworks.com/matlabcentral/fileexchange/6543">File Exchange</a>  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=EVW2A4G2HBVAU">Donate via PayPal</a>

   Markus Buehren
   Last modified 05.07.2011

ASSIGNMENTSUBOPTIMAL1    Compute suboptimal assignment
   ASSIGNMENTSUBOPTIMAL1(DISTMATRIX) computes a suboptimal assignment
   (minimum overall costs) for the given rectangular distance or cost
   matrix, for example the assignment of tracks (in rows) to observations
   (in columns). The result is a column vector containing the assigned
   column number in each row (or 0 if no assignment could be done).

   [ASSIGNMENT, COST] = ASSIGNMENTSUBOPTIMAL1(DISTMATRIX) returns the 
   assignment vector and the overall cost. 

   The algorithm is designed for distance matrices with many forbidden and 
   singly validated assignments (rows or columns containing only one
   finite element). The algorithm first searches the matrix for singly
   validated columns and rejects all assignments with multiply validated
   rows. Afterwards, singly validated rows are searched and assignments to
   multiply validated columns are rejected. Then, for each row that
   validates only with singly validated columns (and the other way
   around), the minimum element is chosen and the assignment is made. If
   there are still assignments open, the minimum element in the distance 
   matrix is searched and the corresponding assignment is made.

   In scenarios without any forbidden assignments, the algorithm reduces
   to the last step, which will provide the same result as ASSIGNMENTOPTIMAL2. 
   If there are only some assignments forbidden, the algorithm will perform
   poorly because singly validated assignments are preferred.

   The last step can still be optimized, see the comments in
   ASSIGNMENTOPTIMAL2.

   <a href="assignment.html">assignment.html</a>  <a href="http://www.mathworks.com/matlabcentral/fileexchange/6543">File Exchange</a>  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=EVW2A4G2HBVAU">Donate via PayPal</a>

   Markus Buehren
   Last modified 05.07.2011

ASSIGNMENTSUBOPTIMAL2    Compute suboptimal assignment
   ASSIGNMENTSUBOPTIMAL2(DISTMATRIX) computes a suboptimal assignment
   (minimum overall costs) for the given rectangular distance or cost
   matrix, for example the assignment of tracks (in rows) to observations
   (in columns). The result is a column vector containing the assigned
   column number in each row (or 0 if no assignment could be done).

   [ASSIGNMENT, COST] = ASSIGNMENTSUBOPTIMAL2(DISTMATRIX) returns the 
   assignment vector and the overall cost.

   The algorithm searches the matrix for the minimum element and makes the
   corresponding row-column assignment. After setting all elements in the
   given row and column to infinity (i.e. forbidden assignment), the
   search procedure is repeated until all assignments are done or only
   infinite values are found.

   This function and the corresponding mex-function can further be
   improved by first sorting all elements instead of searching for the
   minimum of all elements many times.

   <a href="assignment.html">assignment.html</a>  <a href="http://www.mathworks.com/matlabcentral/fileexchange/6543">File Exchange</a>  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=EVW2A4G2HBVAU">Donate via PayPal</a>

   Markus Buehren
   Last modified 05.07.2011
*/

/* MEX-C code references:
	ref: https://www.mathworks.com/help/matlab/matlab_external/c-mex-source-file.html
	ref: https://www.mathworks.com/help/matlab/matlab_external/cpp-mex-api.html
	ref: https://www.mathworks.com/help/matlab/cc-mx-matrix-library.html
	ref: https://stackoverflow.com/questions/3437404/min-and-max-in-c
*/

#define CHECK_FOR_INF
#define ONE_INDEXING

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#if !defined(ABS)
#define ABS(A) ((A > 0) ? (A) : 0)
#endif

/* assignmentoptimal */
void assignmentoptimal(double *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);
void buildassignmentvector(double *assignment, bool *starMatrix, int nOfRows, int nOfColumns);
void computeassignmentcost(double *assignment, double *cost, double *distMatrix, int nOfRows);
void step2a(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
void step2b(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
void step3 (double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
void step4 (double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col);
void step5 (double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);

/* assignmentsuboptimal1 */
void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);

/* assignmentsuboptimal2 */
void assignmentsuboptimal2(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	/* Check arguments */
    if (nrhs < 2 || nrhs > 4){
        mexErrMsgTxt("Two input required and at most four input acceptable");
    }

    if (nlhs > 3){
        mexErrMsgTxt("Maximum three output supported");
    }

    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
    {
        mexErrMsgTxt("Input must be non-complex double elements");
    }

    if (mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1])){
        mexErrMsgTxt("Input must be same size");
    }

    if (MIN(mxGetM(prhs[0]), mxGetN(prhs[0])) != 1 || MIN(mxGetM(prhs[1]), mxGetN(prhs[1])) != 1){
        mexErrMsgTxt("Input must be an array");
    }

	if (nrhs >= 3){
		if (!mxIsScalar(prhs[2]))
			mexErrMsgTxt("3rd input must be a scalar");
		if (nrhs >= 4){
			if (!mxIsScalar(prhs[3]))
				mexErrMsgTxt("4th input must be a scalar");
		}
	}

	/* convert to column  */
	bool flag_row = (mxGetN(prhs[1]) != 1);

    /* Pointers needed */
	double *eps, *eps_hat, cirbnd, method;
	double *eps_paired, *assignment, *cost;

	/* Get data */
    const size_t num_ele    = mxGetNumberOfElements(prhs[0]);

	eps 	= mxGetPr(prhs[0]);
	eps_hat = mxGetPr(prhs[1]);

	if (nrhs >= 3)
		cirbnd = mxGetPr(prhs[2])[0];
	else
		cirbnd = -1;

	if (nrhs >= 4)
		method = mxGetPr(prhs[3])[0];
	else
		method = 0;

	// the following code seems to have some issue when input is double but get int
	// cirbnd = (nrhs >= 3) ? mxGetPr(prhs[2])[0] : -1;
	// method 	= (nrhs >= 4) ? mxGetPr(prhs[3])[0] : 0;

	// check method:
	if (method < -1e-9 || method > 2+1e-9)
		mexErrMsgTxt("method should be 0, 1, 2");

	/* Compute distance matrix */
	double *distMatrix = (double*)mxMalloc(num_ele * num_ele * sizeof(double));

	if (cirbnd > 0){
		for (int i = 0; i < num_ele; i++){
			for (int j = 0; j < num_ele; j++){
				double d = ABS(eps[i] - eps_hat[j]);
				distMatrix[i + num_ele*j] = MIN(d, 2 * cirbnd - d);
			}
		}
	}
	else{
		for (int i = 0; i < num_ele; i++){
			for (int j = 0; j < num_ele; j++){
				distMatrix[i + num_ele*j] = ABS(eps[i] - eps_hat[j]);
			}
		}
	}	
	
	/* Output arguments */
	if (flag_row){
		plhs[0]    = mxCreateDoubleMatrix(1, num_ele, mxREAL);
		plhs[1]    = mxCreateDoubleMatrix(1, num_ele, mxREAL);
	}
	else{
		plhs[0]    = mxCreateDoubleMatrix(num_ele, 1, mxREAL);
		plhs[1]    = mxCreateDoubleMatrix(num_ele, 1, mxREAL);
	}
	plhs[2]    = mxCreateDoubleScalar(0);

    /* Assign pointers to the various parameters */
	eps_paired = mxGetPr(plhs[0]);
	assignment = mxGetPr(plhs[1]);
	cost       = mxGetPr(plhs[2]);
	
	/* Call C-function */
	if (method < 0.5)
		assignmentoptimal(assignment, cost, distMatrix, num_ele, num_ele);
	else if (method < 1.5)
		assignmentsuboptimal1(assignment, cost, distMatrix, num_ele, num_ele);
	else
		assignmentsuboptimal2(assignment, cost, distMatrix, num_ele, num_ele);

    /* Assigen indices to eps_hat */
    for (int i = 0; i < num_ele; i++)
        eps_paired[i] = eps_hat[(mxInt32)assignment[i] - 1];

	/* Delete distance matrix */
	mxFree(distMatrix);
}

/*======================================================*/
/*================== assignmentoptimal =================*/
/*======================================================*/
void assignmentoptimal(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{
	double *distMatrix, *distMatrixTemp, *distMatrixEnd, *columnEnd, value, minValue;
	bool *coveredColumns, *coveredRows, *starMatrix, *newStarMatrix, *primeMatrix;
	int nOfElements, minDim, row, col;
#ifdef CHECK_FOR_INF
	bool infiniteValueFound;
	double maxFiniteValue, infValue;
#endif
	
	/* initialization */
	*cost = 0;
	for(row=0; row<nOfRows; row++)
#ifdef ONE_INDEXING
		assignment[row] =  0.0;
#else
		assignment[row] = -1.0;
#endif
	
	/* generate working copy of distance Matrix */
	/* check if all matrix elements are positive */
	nOfElements   = nOfRows * nOfColumns;
	distMatrix    = (double *)mxMalloc(nOfElements * sizeof(double));
	distMatrixEnd = distMatrix + nOfElements;
	for(row=0; row<nOfElements; row++)
	{
		value = distMatrixIn[row];
		if(mxIsFinite(value) && (value < 0))
			mexErrMsgTxt("All matrix elements have to be non-negative.");
		distMatrix[row] = value;
	}

#ifdef CHECK_FOR_INF
	/* check for infinite values */
	maxFiniteValue     = -1;
	infiniteValueFound = false;
	
	distMatrixTemp = distMatrix;
	while(distMatrixTemp < distMatrixEnd)
	{
		value = *distMatrixTemp++;
		if(mxIsFinite(value))
		{
			if(value > maxFiniteValue)
				maxFiniteValue = value;
		}
		else
			infiniteValueFound = true;
	}
	if(infiniteValueFound)
	{
		if(maxFiniteValue == -1) /* all elements are infinite */
			return;
		
		/* set all infinite elements to big finite value */
		if(maxFiniteValue > 0)
			infValue = 10 * maxFiniteValue * nOfElements;
		else
			infValue = 10;
		distMatrixTemp = distMatrix;
		while(distMatrixTemp < distMatrixEnd)
			if(mxIsInf(*distMatrixTemp++))
				*(distMatrixTemp-1) = infValue;
	}
#endif
				
	/* memory allocation */
	coveredColumns = (bool *)mxCalloc(nOfColumns,  sizeof(bool));
	coveredRows    = (bool *)mxCalloc(nOfRows,     sizeof(bool));
	starMatrix     = (bool *)mxCalloc(nOfElements, sizeof(bool));
	primeMatrix    = (bool *)mxCalloc(nOfElements, sizeof(bool));
	newStarMatrix  = (bool *)mxCalloc(nOfElements, sizeof(bool)); /* used in step4 */

	/* preliminary steps */
	if(nOfRows <= nOfColumns)
	{
		minDim = nOfRows;
		
		for(row=0; row<nOfRows; row++)
		{
			/* find the smallest element in the row */
			distMatrixTemp = distMatrix + row;
			minValue = *distMatrixTemp;
			distMatrixTemp += nOfRows;			
			while(distMatrixTemp < distMatrixEnd)
			{
				value = *distMatrixTemp;
				if(value < minValue)
					minValue = value;
				distMatrixTemp += nOfRows;
			}
			
			/* subtract the smallest element from each element of the row */
			distMatrixTemp = distMatrix + row;
			while(distMatrixTemp < distMatrixEnd)
			{
				*distMatrixTemp -= minValue;
				distMatrixTemp += nOfRows;
			}
		}
		
		/* Steps 1 and 2a */
		for(row=0; row<nOfRows; row++)
			for(col=0; col<nOfColumns; col++)
				if(distMatrix[row + nOfRows*col] == 0)
					if(!coveredColumns[col])
					{
						starMatrix[row + nOfRows*col] = true;
						coveredColumns[col]           = true;
						break;
					}
	}
	else /* if(nOfRows > nOfColumns) */
	{
		minDim = nOfColumns;
		
		for(col=0; col<nOfColumns; col++)
		{
			/* find the smallest element in the column */
			distMatrixTemp = distMatrix     + nOfRows*col;
			columnEnd      = distMatrixTemp + nOfRows;
			
			minValue = *distMatrixTemp++;			
			while(distMatrixTemp < columnEnd)
			{
				value = *distMatrixTemp++;
				if(value < minValue)
					minValue = value;
			}
			
			/* subtract the smallest element from each element of the column */
			distMatrixTemp = distMatrix + nOfRows*col;
			while(distMatrixTemp < columnEnd)
				*distMatrixTemp++ -= minValue;
		}
		
		/* Steps 1 and 2a */
		for(col=0; col<nOfColumns; col++)
			for(row=0; row<nOfRows; row++)
				if(distMatrix[row + nOfRows*col] == 0)
					if(!coveredRows[row])
					{
						starMatrix[row + nOfRows*col] = true;
						coveredColumns[col]           = true;
						coveredRows[row]              = true;
						break;
					}
		for(row=0; row<nOfRows; row++)
			coveredRows[row] = false;
		
	}	
	
	/* move to step 2b */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

	/* compute cost and remove invalid assignments */
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows);
	
	/* free allocated memory */
	mxFree(distMatrix);
	mxFree(coveredColumns);
	mxFree(coveredRows);
	mxFree(starMatrix);
	mxFree(primeMatrix);
	mxFree(newStarMatrix);

	return;
}

/********************************************************/
void buildassignmentvector(double *assignment, bool *starMatrix, int nOfRows, int nOfColumns)
{
	int row, col;
	
	for(row=0; row<nOfRows; row++)
		for(col=0; col<nOfColumns; col++)
			if(starMatrix[row + nOfRows*col])
			{
#ifdef ONE_INDEXING
				assignment[row] = col + 1; /* MATLAB-Indexing */
#else
				assignment[row] = col;
#endif
				break;
			}
}

/********************************************************/
void computeassignmentcost(double *assignment, double *cost, double *distMatrix, int nOfRows)
{
	int row, col;
#ifdef CHECK_FOR_INF
	double value;
#endif
	
	for(row=0; row<nOfRows; row++)
	{
#ifdef ONE_INDEXING
		col = assignment[row]-1; /* MATLAB-Indexing */
#else
		col = assignment[row];
#endif

		if(col >= 0)
		{
#ifdef CHECK_FOR_INF
			value = distMatrix[row + nOfRows*col];
			if(mxIsFinite(value))
				*cost += value;
			else
#ifdef ONE_INDEXING
				assignment[row] =  0.0;
#else
				assignment[row] = -1.0;
#endif

#else
			*cost += distMatrix[row + nOfRows*col];
#endif
		}
	}
}

/********************************************************/
void step2a(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool *starMatrixTemp, *columnEnd;
	int col;
	
	/* cover every column containing a starred zero */
	for(col=0; col<nOfColumns; col++)
	{
		starMatrixTemp = starMatrix     + nOfRows*col;
		columnEnd      = starMatrixTemp + nOfRows;
		while(starMatrixTemp < columnEnd){
			if(*starMatrixTemp++)
			{
				coveredColumns[col] = true;
				break;
			}
		}	
	}

	/* move to step 3 */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void step2b(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	int col, nOfCoveredColumns;
	
	/* count covered columns */
	nOfCoveredColumns = 0;
	for(col=0; col<nOfColumns; col++)
		if(coveredColumns[col])
			nOfCoveredColumns++;
			
	if(nOfCoveredColumns == minDim)
	{
		/* algorithm finished */
		buildassignmentvector(assignment, starMatrix, nOfRows, nOfColumns);
	}
	else
	{
		/* move to step 3 */
		step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}
	
}

/********************************************************/
void step3(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	bool zerosFound;
	int row, col, starCol;

	zerosFound = true;
	while(zerosFound)
	{
		zerosFound = false;		
		for(col=0; col<nOfColumns; col++)
			if(!coveredColumns[col])
				for(row=0; row<nOfRows; row++)
					if((!coveredRows[row]) && (distMatrix[row + nOfRows*col] == 0))
					{
						/* prime zero */
						primeMatrix[row + nOfRows*col] = true;
						
						/* find starred zero in current row */
						for(starCol=0; starCol<nOfColumns; starCol++)
							if(starMatrix[row + nOfRows*starCol])
								break;
						
						if(starCol == nOfColumns) /* no starred zero found */
						{
							/* move to step 4 */
							step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						}
						else
						{
							coveredRows[row]        = true;
							coveredColumns[starCol] = false;
							zerosFound              = true;
							break;
						}
					}
	}
	
	/* move to step 5 */
	step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void step4(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col)
{	
	int n, starRow, starCol, primeRow, primeCol;
	int nOfElements = nOfRows*nOfColumns;
	
	/* generate temporary copy of starMatrix */
	for(n=0; n<nOfElements; n++)
		newStarMatrix[n] = starMatrix[n];
	
	/* star current zero */
	newStarMatrix[row + nOfRows*col] = true;

	/* find starred zero in current column */
	starCol = col;
	for(starRow=0; starRow<nOfRows; starRow++)
		if(starMatrix[starRow + nOfRows*starCol])
			break;

	while(starRow<nOfRows)
	{
		/* unstar the starred zero */
		newStarMatrix[starRow + nOfRows*starCol] = false;
	
		/* find primed zero in current row */
		primeRow = starRow;
		for(primeCol=0; primeCol<nOfColumns; primeCol++)
			if(primeMatrix[primeRow + nOfRows*primeCol])
				break;
								
		/* star the primed zero */
		newStarMatrix[primeRow + nOfRows*primeCol] = true;
	
		/* find starred zero in current column */
		starCol = primeCol;
		for(starRow=0; starRow<nOfRows; starRow++)
			if(starMatrix[starRow + nOfRows*starCol])
				break;
	}	

	/* use temporary copy as new starMatrix */
	/* delete all primes, uncover all rows */
	for(n=0; n<nOfElements; n++)
	{
		primeMatrix[n] = false;
		starMatrix[n]  = newStarMatrix[n];
	}
	for(n=0; n<nOfRows; n++)
		coveredRows[n] = false;
	
	/* move to step 2a */
	step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void step5(double *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim)
{
	double h, value;
	int row, col;
	
	/* find smallest uncovered element h */
	h = mxGetInf();	
	for(row=0; row<nOfRows; row++)
		if(!coveredRows[row])
			for(col=0; col<nOfColumns; col++)
				if(!coveredColumns[col])
				{
					value = distMatrix[row + nOfRows*col];
					if(value < h)
						h = value;
				}
	
	/* add h to each covered row */
	for(row=0; row<nOfRows; row++)
		if(coveredRows[row])
			for(col=0; col<nOfColumns; col++)
				distMatrix[row + nOfRows*col] += h;
	
	/* subtract h from each uncovered column */
	for(col=0; col<nOfColumns; col++)
		if(!coveredColumns[col])
			for(row=0; row<nOfRows; row++)
				distMatrix[row + nOfRows*col] -= h;
	
	/* move to step 3 */
	step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/*======================================================*/
/*================ assignmentsuboptimal1 ===============*/
/*======================================================*/
void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{
	bool infiniteValueFound, finiteValueFound, repeatSteps, allSinglyValidated, singleValidationFound;
	int n, row, col, tmpRow, tmpCol, nOfElements;
	int *nOfValidObservations, *nOfValidTracks;
	double value, minValue, *distMatrix, inf;
	
	inf = mxGetInf();
	
	/* make working copy of distance Matrix */
	nOfElements   = nOfRows * nOfColumns;
	distMatrix    = (double *)mxMalloc(nOfElements * sizeof(double));
	for(n=0; n<nOfElements; n++)
		distMatrix[n] = distMatrixIn[n];
	
	/* initialization */
	*cost = 0;
#ifdef ONE_INDEXING
	for(row=0; row<nOfRows; row++)
		assignment[row] =  0.0;
#else
	for(row=0; row<nOfRows; row++)
		assignment[row] = -1.0;
#endif
	
	/* allocate memory */
	nOfValidObservations  = (int *)mxCalloc(nOfRows,    sizeof(int));
	nOfValidTracks        = (int *)mxCalloc(nOfColumns, sizeof(int));
		
	/* compute number of validations */
	infiniteValueFound = false;
	finiteValueFound  = false;
	for(row=0; row<nOfRows; row++)
		for(col=0; col<nOfColumns; col++)
			if(mxIsFinite(distMatrix[row + nOfRows*col]))
			{
				nOfValidTracks[col]       += 1;
				nOfValidObservations[row] += 1;
				finiteValueFound = true;
			}
			else
				infiniteValueFound = true;
				
	if(infiniteValueFound)
	{
		if(!finiteValueFound)
			return;
			
		repeatSteps = true;
		
		while(repeatSteps)
		{
			repeatSteps = false;

			/* step 1: reject assignments of multiply validated tracks to singly validated observations		 */
			for(col=0; col<nOfColumns; col++)
			{
				singleValidationFound = false;
				for(row=0; row<nOfRows; row++)
					if(mxIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidObservations[row] == 1))
					{
						singleValidationFound = true;
						break;
					}
					
				if(singleValidationFound)
				{
					for(row=0; row<nOfRows; row++)
						if((nOfValidObservations[row] > 1) && mxIsFinite(distMatrix[row + nOfRows*col]))
						{
							distMatrix[row + nOfRows*col] = inf;
							nOfValidObservations[row] -= 1;							
							nOfValidTracks[col]       -= 1;	
							repeatSteps = true;				
						}
					}
			}
			
			/* step 2: reject assignments of multiply validated observations to singly validated tracks */
			if(nOfColumns > 1)			
			{	
				for(row=0; row<nOfRows; row++)
				{
					singleValidationFound = false;
					for(col=0; col<nOfColumns; col++)
						if(mxIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidTracks[col] == 1))
						{
							singleValidationFound = true;
							break;
						}
						
					if(singleValidationFound)
					{
						for(col=0; col<nOfColumns; col++)
							if((nOfValidTracks[col] > 1) && mxIsFinite(distMatrix[row + nOfRows*col]))
							{
								distMatrix[row + nOfRows*col] = inf;
								nOfValidObservations[row] -= 1;
								nOfValidTracks[col]       -= 1;
								repeatSteps = true;								
							}
						}
				}
			}
		} /* while(repeatSteps) */
	
		/* for each multiply validated track that validates only with singly validated  */
		/* observations, choose the observation with minimum distance */
		for(row=0; row<nOfRows; row++)
		{
			if(nOfValidObservations[row] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(col=0; col<nOfColumns; col++)
				{
					value = distMatrix[row + nOfRows*col];
					if(mxIsFinite(value))
					{
						if(nOfValidTracks[col] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidTracks[col] == 1) && (value < minValue))
						{
							tmpCol   = col;
							minValue = value;
						}
					}
				}
				
				if(allSinglyValidated)
				{
	#ifdef ONE_INDEXING
					assignment[row] = tmpCol + 1;
	#else
					assignment[row] = tmpCol;
	#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*tmpCol] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[row + nOfRows*n] = inf;
				}
			}
		}

		/* for each multiply validated observation that validates only with singly validated  */
		/* track, choose the track with minimum distance */
		for(col=0; col<nOfColumns; col++)
		{
			if(nOfValidTracks[col] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(row=0; row<nOfRows; row++)
				{
					value = distMatrix[row + nOfRows*col];
					if(mxIsFinite(value))
					{
						if(nOfValidObservations[row] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidObservations[row] == 1) && (value < minValue))
						{
							tmpRow   = row;
							minValue = value;
						}
					}
				}
				
				if(allSinglyValidated)
				{
	#ifdef ONE_INDEXING
					assignment[tmpRow] = col + 1;
	#else
					assignment[tmpRow] = col;
	#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*col] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[tmpRow + nOfRows*n] = inf;
				}
			}
		}	
	} /* if(infiniteValueFound) */
	
	
	/* now, recursively search for the minimum element and do the assignment */
	while(true)
	{
		/* find minimum distance observation-to-track pair */
		minValue = inf;
		for(row=0; row<nOfRows; row++)
			for(col=0; col<nOfColumns; col++)
			{
				value = distMatrix[row + nOfRows*col];
				if(mxIsFinite(value) && (value < minValue))
				{
					minValue = value;
					tmpRow   = row;
					tmpCol   = col;
				}
			}
		
		if(mxIsFinite(minValue))
		{
#ifdef ONE_INDEXING
			assignment[tmpRow] = tmpCol+ 1;
#else
			assignment[tmpRow] = tmpCol;
#endif
			*cost += minValue;
			for(n=0; n<nOfRows; n++)
				distMatrix[n + nOfRows*tmpCol] = inf;
			for(n=0; n<nOfColumns; n++)
				distMatrix[tmpRow + nOfRows*n] = inf;			
		}
		else
			break;
			
	} /* while(true) */
	
	/* free allocated memory */
	mxFree(nOfValidObservations);
	mxFree(nOfValidTracks);


}

/*======================================================*/
/*================ assignmentsuboptimal2 ===============*/
/*======================================================*/
void assignmentsuboptimal2(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{
	int n, row, col, tmpRow, tmpCol, nOfElements;
	double value, minValue, *distMatrix, inf;

	inf = mxGetInf();
	
	/* make working copy of distance Matrix */
	nOfElements   = nOfRows * nOfColumns;
	distMatrix    = (double *)mxMalloc(nOfElements * sizeof(double));
	for(n=0; n<nOfElements; n++)
		distMatrix[n] = distMatrixIn[n];
	
	/* initialization */
	*cost = 0;
	for(row=0; row<nOfRows; row++)
#ifdef ONE_INDEXING
		assignment[row] =  0.0;
#else
		assignment[row] = -1.0;
#endif
	
	/* recursively search for the minimum element and do the assignment */
	while(true)
	{
		/* find minimum distance observation-to-track pair */
		minValue = inf;
		for(row=0; row<nOfRows; row++)
			for(col=0; col<nOfColumns; col++)
			{
				value = distMatrix[row + nOfRows*col];
				if(mxIsFinite(value) && (value < minValue))
				{
					minValue = value;
					tmpRow   = row;
					tmpCol   = col;
				}
			}
		
		if(mxIsFinite(minValue))
		{
#ifdef ONE_INDEXING
			assignment[tmpRow] = tmpCol+ 1;
#else
			assignment[tmpRow] = tmpCol;
#endif
			*cost += minValue;
			for(n=0; n<nOfRows; n++)
				distMatrix[n + nOfRows*tmpCol] = inf;
			for(n=0; n<nOfColumns; n++)
				distMatrix[tmpRow + nOfRows*n] = inf;			
		}
		else
			break;
			
	} /* while(true) */
	
	mxFree(distMatrix);
}
