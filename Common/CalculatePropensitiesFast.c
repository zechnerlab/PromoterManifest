#include "mex.h"
#include "math.h"


int Factorial(int x)
{
    if (x == 0)
    {
        return 1;
    }
    
    return x * Factorial(x - 1);
}

int NChooseK(int n, int k)
{
    /* return Factorial(n) / (Factorial(k) * Factorial(n - k)) */;
    
    int i, prod;
    
    prod = 1;
    
    for (i=1; i<=k; i++)
    {
        prod = prod * (n - k + i) / i;
    }
    
    return prod;
}


void CalculatePropensities(double *h, double *g, double *x,
        double *c, double *pre, mwSize numSpecies, mwSize numReactions,
        mwSize numStates)
{
    
    mwSize i;
    mwSize k;
    mwSize l;
    mwSize propIdx;
    
    int stochProduct;
    int preVal;
    int stateValue;

    for (l=0; l<numStates; l++)
    {
        for (i=0; i<numReactions; i++)
        {
            stochProduct = 1;
            
            for (k=0; k<numSpecies; k++)
            {
                preVal = (int)pre[numReactions * k + i];
                
                stateValue = (int)x[l*numSpecies + k];
                
                if (stateValue >= preVal)
                {
                    stochProduct = stochProduct * NChooseK(stateValue, preVal);
                }
                else
                {
                    stochProduct = 0;
                    break;
                }
            }
            
            propIdx = l*numReactions + i;
            g[propIdx] = stochProduct;
            h[propIdx] = c[i] * stochProduct;
        }
    }
    
    
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *h, *g, *c;
  double *x, *pre;
  mwSize numSpecies, numReactions, numStates;
  
  /* Check for proper number of arguments. */
  if(nrhs!=3) {
    mexErrMsgTxt("Three inputs required.");
  } else if(nlhs>2) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  /* The input must be a noncomplex scalar double.*/
  numSpecies = mxGetN(prhs[2]);
  numReactions = mxGetM(prhs[2]);
  numStates = mxGetN(prhs[0]);
/*   if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !(mrows==1 && ncols==1) ) {
        mexErrMsgTxt("Input must be a noncomplex scalar double.");
     } */
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(numReactions,numStates, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(numReactions,numStates, mxREAL);
  
  /* Assign pointers to each input and output. */
  h = mxGetPr(plhs[0]);
  g = mxGetPr(plhs[1]);
  
  x = mxGetPr(prhs[0]);
  c = mxGetPr(prhs[1]);
  pre = mxGetPr(prhs[2]);
  
  /* Call the timestwo subroutine. */
  CalculatePropensities(h, g, x, c, pre, numSpecies, numReactions, numStates);
}
