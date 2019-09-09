#include "mex.h"
#include "math.h"
#include "string.h"

void display(const double* in);
void copy(double* res, const double* vec, const int n);
void addto(double* res, const double* vec, const int n);
void add(double* res, const double* vec1, const double* vec2, const int n);
void zero(double* res, const int n);
void scale(double* out, double* in, const double factor, const int n);
void transpose(double* res, const double* mat, const int n, const int m);
void multMatTMat(double* res, const double* mat1, const double* mat2,
				   const int r1, const int c1, const int c2);
void multMatMat(double* res, const double* mat1, const double* mat2,
				   const int r1, const int c1, const int c2);
void multVecTens(double* res, const double* vec, const double* tens,
				   const int r, const int c, const int p);
void multVecSymTens(double* res, const double* vec, const double* tens,
				   const int r, const int c, const int p);
void multVecDiagTens(double* res, const double* vec, const double* tens,
				   const int n);
void multMatSymMat(double* res, const double* mat1, const double* mat2,
				   const int r1, const int c1, const int c2);
int mju_solverChol(double* res, double* mat, const double* vec, 
				   const int n, const int nvec, const char flg_factorize);

/* [div, Vx, Vxx, l, L, dV] = back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lam,reg) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* ======== Declarations */
   double *div, *l, *L, *Vx, *Vxx, *dV,
           *cx, *cu, *cxx, *cxu, *cuu, *fx, *fu, *fxx, *fxu, *fuu, *lam, *reg,
           *Qu, *Qux, *Quu, *QuxF, *QuuF, *temp1, *temp2,
           *l_k, *L_k, *Vx_k, *Vxx_k, *cx_k, *cu_k, *cxx_k, *cxu_k, *cuu_k, 
           *fx_k, *fu_k, *fxx_k, *fxu_k, *fuu_k, 
           *Vx_p, *Vxx_p;
   
   mwSize n, N, m, get_fxx;

   mwSize *dims;
   
   const mwSize *dims1;
   
   //char *chN = "N", *chT = "T", *chU = "U", *chL = "L", *chR = "R";
   
   int i,j,k,nn,nm;
   
   /* ======== Dimension checks, initialization, memory allocation */
   if ((nrhs != 12)){
      mexErrMsgTxt("The number of input arguments must be 12.");
   }
   
   dims     = (mwSize *) mxCalloc(3 , sizeof( mwSize ));
   
   /* RHS and related */
   if (mxGetNumberOfDimensions(prhs[0]) == 2){
      n    = mxGetM(prhs[0]);
      N    = mxGetN(prhs[0]);
      m    = mxGetM(prhs[1]);  
   }else{
      dims1 = mxGetDimensions(prhs[0]);
      n    = dims1[1];
      N    = dims1[2];
      dims1 = mxGetDimensions(prhs[1]);
      m    = dims1[1];
   }
   
   cx   = mxGetPr(prhs[0]);
   cu   = mxGetPr(prhs[1]);
   cxx  = mxGetPr(prhs[2]);
   cxu  = mxGetPr(prhs[3]);
   cuu  = mxGetPr(prhs[4]);
   fx   = mxGetPr(prhs[5]);
   fu   = mxGetPr(prhs[6]);
   fxx  = mxGetPr(prhs[7]);
   get_fxx = mxGetNumberOfDimensions(prhs[7]);   
   fxu  = mxGetPr(prhs[8]);
   fuu  = mxGetPr(prhs[9]);   
   lam  = mxGetPr(prhs[10]);
   reg  = mxGetPr(prhs[11]);
   
   /* LHS  and related */   
   plhs[0]  = mxCreateDoubleScalar(0); /*allocate for div */
   div      = mxGetPr(plhs[0]);
   
   dims[0]  = n;
   dims[1]  = N;
   plhs[1]  = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); /*allocate for Vx */
   Vx       = mxGetPr(plhs[1]);

   dims[0]  = n;
   dims[1]  = n;
   dims[2]  = N;
   plhs[2]     = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); /*allocate for Vxx */
   Vxx         = mxGetPr(plhs[2]);

   dims[0]  = m;
   dims[1]  = N-1;
   plhs[3]  = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); /*allocate for l */   
   l        = mxGetPr(plhs[3]);

   dims[0]  = m;
   dims[1]  = n;
   dims[2]  = N-1;
   plhs[4]  = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); /*allocate for L */   
   L        = mxGetPr(plhs[4]);

   dims[0]  = 1;
   dims[1]  = 2;
   plhs[5]  = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); /*allocate for dV */    
   dV       = mxGetPr(plhs[5]);
   
   /* allocate for workspace variables */
   Qu          = mxCalloc(m,     sizeof( double ));
   Qux         = mxCalloc(m*n,   sizeof( double ));
   Quu         = mxCalloc(m*m,   sizeof( double ));
   QuxF        = mxCalloc(m*n,   sizeof( double ));
   QuuF        = mxCalloc(m*m,   sizeof( double ));
   temp1       = mxCalloc(2*n*n, sizeof( double ));
   temp2       = mxCalloc(2*n*n, sizeof( double ));   
   
   /* set final Vx and Vxx to final costs */
   for (i=0; i<n; i++){
      Vx[n*(N-1) + i] = cx[n*(N-1) + i];
      for (j=0; j<n; j++)
         Vxx[n*n*(N-1) + i*n + j] = cxx[n*n*(N-1) + i*n + j]; 
   }
   
   /* main loop */
   for (k=N-2; k>-1; k--){
      
// indexing shortcuts for readabillity
      cx_k   = cx  + n*k;
      cu_k   = cu  + m*k;
      cxx_k  = cxx + n*n*k;
      cxu_k  = cxu + n*m*k;
      cuu_k  = cuu + m*m*k;
      fx_k   = fx  + n*n*k;
      fu_k   = fu  + n*m*k;
      fxu_k  = fxu + n*n*m*k;
      fuu_k  = fuu + n*m*m*k;      
      l_k    = l   + m*k;
      L_k    = L   + m*n*k;
      Vx_k   = Vx  + n*k;
      Vxx_k  = Vxx + n*n*k;
      Vx_p   = Vx  + n*(k+1);
      Vxx_p  = Vxx + n*n*(k+1); 
      
      if (get_fxx == 4)
         fxx_k  = fxx + n*n*n*k;
      else if (get_fxx == 3)
         fxx_k  = fxx + n*n*k;      
      
// compute approximation to Q(x,u)      
      
      //    Qu  = cu(:,k) + fu(:,:,k)'*Vx(:,k+1); 
      multMatTMat(Qu, fu_k, Vx_p, n, m, 1);
      addto(Qu, cu_k, m);
      
      //    Qux = cxu(:,:,k)'  + fu(:,:,k)'*Vxx(:,:,k+1)*fx(:,:,k) + p13mm(Vx(:,k+1),fxu(:,:,:,k));
      multMatTMat(temp1, fu_k, Vxx_p, n, m, n);
      multMatMat(Qux, temp1, fx_k, m, n, n);
      transpose(temp2,cxu_k,n,m);
      addto(Qux, temp2, m*n);
      multVecTens(temp2, Vx_p, fxu_k, n, n, m);
      addto(Qux, temp2, m*n);
      
      //    Quu = cuu(:,:,k)   + fu(:,:,k)'*Vxx(:,:,k+1)*fu(:,:,k) + p13mm(Vx(:,k+1),fuu(:,:,:,k));
      multMatMat(Quu, temp1, fu_k, m, n, m);
      addto(Quu, cuu_k, m*m);
      multVecSymTens(temp1, Vx_p, fuu_k, n, m, m);
      addto(Quu, temp1, m*m);
      
      //    VxxF = (Vxx(:,:,k+1) + lambda*eye(n)*(regType == 2));
      //    QuxF = cxu(:,:,k)' + fu(:,:,k)'*VxxF*fx(:,:,k) + p13mm(Vx(:,k+1),fxu(:,:,:,k));         
      //    QuuF = cuu(:,:,k)  + fu(:,:,k)'*VxxF*fu(:,:,k) + lambda*eye(m)*(regType == 1);      
      copy(QuxF, Qux, n*m);
      copy(QuuF, Quu, m*m);
      if ((int)*reg == 2){
         multMatTMat(temp1, fu_k, fx_k, n, m, n);
         scale(temp1, temp1, *lam, m*n);
         addto(QuxF, temp1, m*n);
         multMatTMat(temp1, fu_k, fu_k, n, m, m);
         scale(temp1, temp1, *lam, m*m);
         addto(QuuF, temp1, m*m);         
      }else
         for( i=0; i<m; i++ )
            QuuF[i*m+i] += *lam;
         
// find control law
      copy(temp1,   Qu,   m);
      copy(temp1+m, QuxF, m*n);
      i = mju_solverChol(temp2, QuuF, temp1, m, n+1, 1);
      if (i!=m){
         *div   = (double)(k+1);
         return;
      }
      scale(l_k, temp2,   -1, m);
      scale(L_k, temp2+m, -1, m*n);
      
// update cost-to-go approximation
      
      // dV       = dV + [l(:,k)'*Qu  .5*l(:,k)'*Quu*l(:,k)];
      multMatTMat(temp1, l_k, Qu, m, 1, 1);
      dV[0] +=  *temp1;
      multMatMat(temp1, Quu, l_k, m, m, 1);
      multMatTMat(temp2, l_k, temp1, m, 1, 1);
      dV[1] +=  .5*(*temp2);
      
      // Vx(:,k) = cx(:,k) + fx(:,:,k)'*Vx(:,k+1) 
      multMatTMat(Vx_k, fx_k, Vx_p, n, n, 1);
      addto(Vx_k, cx_k, n);
      // ... + L(:,:,k)'*(Quu*l(:,k) + Qu)
      add(temp2, temp1, Qu, m); //reuse temp1
      multMatTMat(temp1, L_k, temp2, m, n, 1);
      addto(Vx_k, temp1, n);
      // ...  + Qux'*l(:,k);
      multMatTMat(temp1, Qux, l_k, m, n, 1);
      addto(Vx_k, temp1, n);
      
      // Vxx(:,:,k) = cxx(:,:,k) + fx(:,:,k)'*Vxx(:,:,k+1)*fx(:,:,k)  
      multMatMat(temp1, Vxx_p, fx_k, n, n, n);
      multMatSymMat(Vxx_k, fx_k, temp1, n, n, n);
      addto(Vxx_k, cxx_k, n*n);
      // ... + L(:,:,k)'*Quu*L(:,:,k)
      multMatMat(temp1, Quu, L_k, m, m, n);
      multMatSymMat(temp2, L_k, temp1, m, n, n);
      addto(Vxx_k, temp2, n*n);
      // .... + L(:,:,k)'*Qux + Qux'*L(:,:,k)
      multMatTMat(temp2, L_k, Qux, m, n, n);
      transpose(temp1, temp2, n, n);
      addto(Vxx_k, temp2, n*n);
      addto(Vxx_k, temp1, n*n);
      // ... + p13mm(Vx(:,k+1),fxx(:,:,:,k));
      if (get_fxx == 4){
         multVecSymTens(temp1, Vx_p, fxx_k, n, n, n);
         addto(Vxx_k, temp1, n*n);
      }else if (get_fxx == 3){
         multVecDiagTens(temp1, Vx_p, fxx_k, n);
         addto(Vxx_k, temp1, n*n);
      }      
      // symmetrize
      transpose(temp1, Vxx_k, n, n);
      addto(Vxx_k, temp1, n*n);
      scale(Vxx_k, Vxx_k, .5, n*n);
   }
   
   /* free workspace variables */
   mxFree(dims);
   mxFree(Qu);
   mxFree(Qux);
   mxFree(Quu);
   mxFree(QuxF);
   mxFree(QuuF);
   mxFree(temp1);
   mxFree(temp2);   
}


void transpose(double* res, const double* mat, const int n, const int m)
{
	int i, j;
   for( i=0; i<m; i++ )
      for( j=0; j<n; j++ )
         res[j*m + i] = mat[i*n + j];
}

void copy(double* res, const double* vec, const int n)
{
	int i;
   for( i=0; i<n; i++ )
      res[i] = vec[i];
}

void addto(double* res, const double* vec, const int n)
{
	int i;
   for( i=0; i<n; i++ )
      res[i] += vec[i];
}

void add(double* res, const double* vec1, const double* vec2, const int n)
{
	int i;
   for( i=0; i<n; i++ )
      res[i] = vec1[i] + vec2[i];
}

void scale(double* out, double* in, const double factor, const int n)
{
	int i;
   for( i=0; i<n; i++ )
      out[i] = factor*in[i];
}

void zero(double* res, const int n)
{
	int i;
   for( i=0; i<n; i++ )
      res[i] = 0;
}

double dot(const double* vec1, const double* vec2, const int n)
{
	int i;
	double res = 0;
	
   for( i=0; i<n; i++ )
      res += vec1[i] * vec2[i];
   
	return res;
}

void multMatSymMat(double* res, const double* mat1, const double* mat2,
				   const int r1, const int c1, const int c2)
{
	int i, j, k;
	double tmp;
	const double *p1, *p2;
   
	for( i=0; i<c2; i++ )
		for( j=i; j<c1; j++ )
		{
			tmp = 0;
			p1 = mat1 + j*r1;
			p2 = mat2 + i*r1;

			for( k=0; k<r1; k++ )
				tmp += p1[k] * p2[k];
			res[i*c1+j] = tmp;
		}
   
	for( i=0; i<c2; i++ ) //Symmetrize
		for( j=0; j<i; j++ )
         res[i*c1+j] = res[j*c1+i];
}

void multMatTMat(double* res, const double* mat1, const double* mat2,
				   const int r1, const int c1, const int c2)
{
	int i, j, k;
	double tmp;
	const double *p1, *p2;
   
	for( i=0; i<c2; i++ )
		for( j=0; j<c1; j++ )
		{
			tmp = 0;
			p1 = mat1 + j*r1;
			p2 = mat2 + i*r1;

			for( k=0; k<r1; k++ )
				tmp += p1[k] * p2[k];
			res[i*c1+j] = tmp;
		}
}


void multMatMat(double* res, const double* mat1, const double* mat2,
				   const int r1, const int c1, const int c2)
{
	int i, j, k;
	double tmp;
	const double *p1, *p2;
   
	for( i=0; i<c2; i++ )
		for( j=0; j<r1; j++ )
		{
			tmp = 0;
			p1 = mat1 + j;
			p2 = mat2 + i*c1;

			for( k=0; k<c1; k++ )
			{
				tmp += (*p1) * (*p2++);
				p1  += r1;
			}
			res[i*r1+j] = tmp;
		}
}


void multVecTens(double* res, const double* vec, const double* tens,
				   const int r, const int c, const int p)
{
	int i, j, k;

	for( i=0; i<c*p; i++ )
      res[i] = 0;   
   
	for( i=0; i<r; i++ )
		for( j=0; j<c; j++ )
			for( k=0; k<p; k++ )
            res[j*p+k] += vec[i] * tens[i + r*c*k + r*j];
}

void multVecSymTens(double* res, const double* vec, const double* tens,
				   const int r, const int c, const int p)
{
	int i, j, k;

	for( i=0; i<c*p; i++ )
      res[i] = 0;   
   
	for( i=0; i<r; i++ )
		for( j=0; j<c; j++ )
			for( k=j; k<p; k++ )
            res[j*p+k] += vec[i] * tens[i + r*c*k + r*j];
   
   for( j=0; j<c; j++ )
      for( k=0; k<j; k++ )
         res[j*p+k] = res[k*p+j];
}

void multVecDiagTens(double* res, const double* vec, const double* tens,
				   const int n)
{
	int i, j, k;

	for( i=0; i<n*n; i++ )
      res[i] = 0;   

	for( i=0; i<n; i++ )
		for( j=0; j<n; j++ )
         res[i*n+i] += vec[j] * tens[n*i + j];   
}

int mju_solverChol(double* res, double* mat, const double* vec, 
				   const int n, const int nvec, const char flg_factorize)
{
	int i, j, ivec, rank=0;
	double tmp;

	// matrix may have been factorized already
	if( flg_factorize )
	{
		// in-place Cholesky factorization, use 'res' to store 1/L(j,j)
		for( j=0; j<n; j++ )
		{
			tmp = mat[j*n+j];
			if( j )
				tmp -= dot(mat+j*n, mat+j*n, j);

			// handle diagonal values below threshold... 
			if( tmp < 0.000000001 )
				res[j] = 0;
			else
			{
				res[j] = (double)(1.0/sqrt(tmp));
				rank++;
			}

			// process off-diagonal entries, modify 'mat'
			for( i=j+1; i<n; i++ )
			{
				mat[i*n+j] -= dot(mat+i*n, mat+j*n, j);
				mat[i*n+j] *= res[j];
			}
		}

		// copy 'res' to diagonal of mat
		for( j=0; j<n; j++ )
			mat[j*n+j] = res[j];
	}


	// solve for multiple vectors via backsubstitution
	for( ivec=0; ivec<nvec; ivec++ )
	{
		// forward substitution: solve  L*res = vec
		copy(res+n*ivec, vec+n*ivec, n);
		for( i=0; i<n; i++ )
		{
			if( i )
				res[i+n*ivec] -= dot(mat+i*n, res+n*ivec, i);

			res[i+n*ivec] *= mat[i*n+i];
		}

		// backward substitution: solve  L'*res = res
		for( i=n-1; i>=0; i-- )
		{
			if( i < n-1 )
				for( j=i+1; j<n; j++ )
					res[i+n*ivec] -= mat[j*n+i] * res[j+n*ivec];

			res[i+n*ivec] *= mat[i*n+i];
		}
	}

	return rank;
}
