/*//////////////////////////////////////////////////////////////////////////
// File name: GravAccelPartialsExterior_mex.c
//
// Author: Siamak Hesar             3/22/2015
//
// This function is adapted from the original code written by Yu Takahashi 
// called AccelInteriorPotential_mex.c. However I implemented Cunningham's approach
// in computing first and second partial derivatives of the potential w.r.t
// the Cartesian coordinates. The original code used Hotine's method.
// 
///////////////////////////////////////////////////////////////////////////
//
// Description:
//
//  This function computes the acceleration, and dynamics 
//  matrix for the exterior potential. It returns the full partial matrix 
//  of acceleration w.r.t state and gravity harmonics partials.
//
// Inputs:
//
//     n_degree   [n.d.]     : Degree of the spherical harmonics
//
//     R_ref      [km]       : Reference distance. usually the reference radius
//
//     mu      	  [km^3/sec^2] : Gravitational Parameter.
//
//     r_vec      [km]       : Field point (spacecraft) vector = [x_sat, y_sat, z_sat]
//
//     Cbar       [n.d.]     : Normalized C spherical harmonics
//
//     Sbar       [n.d.]     : Normalized S spherical harmonics
//
// Outputs:
//
//     Accel              : Acceleration by basis \bar{K}_{nm}^i (normalized)
//
//     Jacobian               : STM by \bar{b}_{nm}^i
//
//     STM_C_bar_1              : STM of Cbar spherical harmonics by \bar{b}_{nm}^i
//
//     STM_S_bar_1              : STM of Sbar spherical harmonics by \bar{b}_{nm}^i
//
// Assumptions/References:
//	- 1: S. V. Bettadpur, "Hotine's geopotential formulation: revisited", Bulletin Geodesique (1995) 69:i35-142
//  - 2: R. A. Werner, "Evaluating Descent and Ascent Trajectories Near Non-Spherical Bodies", Technical Support Package
//	- 3: L. E. Cunningham, "On the computation of the spherical harmonic terms needed during the numerical integration of the orbital motion of an artificial satellite"
//	
//  - Mathematical Formulation
//
//    ~ Note the following definitions
//
//       (1) b_{n,m}^e       = (R_ref/r)^{n+1} * Pnm * [cos(m*lambda); sin(m*lambda)]
//
//       (2) \bar{b}_{n,m}^e = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!} * (R_ref/r)^{n+1} * Pnm * [cos(m*lambda); sin(m*lambda)]
//
//       (3) c_{n,m}^e       = ( 2 - \delta_{0,m} )*(n - m)!/(n + m)!*(r'/R_ref)^n * Pnm * [cos(m*lambda'); sin(m*lambda')]
//
//       (4) \bar{c}_{n,m}^e = \sqrt{ (2 - \delta_{0,m})*(n - m)!/( (2n + 1)*(n + m)! ) } *(r'/R_ref)^n * Pnm * [cos(m*lambda'); sin(m*lambda')]
//
//      where ' indicates the parameters of the differential mass.
//
//    ~ Note that these expressions can be considered as imaginary numbers. That is,
//
//       (1) b_{n,m}^i       = (R_ref/r)^{n+1} * Pnm * e^{i*m*lambda}
//
//       (2) \bar{b}_{n,m}^i = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!} * (R_ref/r)^{n+1} * Pnm * e^{i*m*lambda}
//
//       (3) c_{n,m}^i       = ( 2 - \delta_{0,m} )*(n - m)!/(n + m)! * (r'/R_ref)^n * Pnm * e^{i*m*lambda'}
//
//       (4) \bar{c}_{n,m}^i = \sqrt{ (2 - \delta_{0,m})*(n - m)!/( (2n + 1)*(n + m)! ) } * (r'/R_ref)^n * Pnm * e^{i*m*lambda'}
//
//    ~ The normalization factor is defined as N, where
//
//        \bar{P}_{n,m} = N*P_{n,m}
//
//      Thus, N is immediately recognized as
//
//         N        = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!}
//
//    ~ The external potetial of the field point is computed as
//
//        (1) U^e       = \frac{G*M_ref*}{R_ref} * \sum^{\infty}_{n = 0} \sum^n_{m = 0} b_{n,m}^i * (1/M_ref) * \int_M c_{n,m}^i dm'
//
//        (2) U^e       = \frac{G*M_ref*}{R_ref} * \sum^{\infty}_{n = 0} \sum^n_{m = 0} \bar{b}_{n,m}^i * (1/M_ref) * \int_M \bar{c}_{n,m}^i dm'
//
//    ~ This function has the following recursive formulae:
//
//      -------------------------------------------------------------------
//
//      (1) Basis function          : b_{0,0}^e   = (R_ref/r) * [1; 0]
//
//      (2) Diagonal recurrences    : b_{n,n}^e   = (2*n - 1) * (R_ref/r) * [x/r, -y/r; y/r, x/r]*b_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : b_{n,n-1}^e = (2*n - 1) * (R_ref/r) * (z/r) * b_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : b_{n,m}^e   = \frac{2*n - 1}{n - m}* (R_ref/r) * (z/r) * b_{n-1,m}^e - \frac{n + m - 1}{n - m}*(R_ref/r)^2*b_{n-2,m}^e
//      -------------------------------------------------------------------
//
//      (1) Basis function          : \bar{b}_{0,0}^e   = (R_ref/r) * [1; 0]
//
//      (2) Diagonal recurrences    : \bar{b}_{n,n}^e   = \sqrt{(1 + \delta_{1,n})*(2n + 1)/(2n)}* (R_ref/r) * [x/r, -y/r; y/r, x/r]*\bar{b}_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : \bar{b}_{n,n-1}^e = \sqrt{2*n - 1}* (R_ref/r) * (z/r)*\bar{b}_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : \bar{b}_{n,m}^e   = \frac{4n^2 - 1}{n^2 - m^2}* (R_ref/r) * (z/r)*\bar{b}_{n-1,m}^e - \frac{(2n + 1)*( (n - 1)^2 - m^2 )}{(2n - 3)*(n^2 - m^2)}*(R_ref/r)^2*\bar{b}_{n-2,m}^e
//
//      -------------------------------------------------------------------
//
//      (1) Basis function          : c_{0,0}^e   = [1; 0]
//
//      (2) Diagonal recurrences    : c_{n,n}^e   = (1 + \delta_{1,n})/(2n)*[x'/R_ref, -y'/R_ref; y'/R_ref, x'/R_ref]*c_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : c_{n,n-1}^e = (z'/R_ref)*c_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : c_{n,m}^e   = \frac{2n - 1}{n + m}*(z'/R_ref)*c_{n-1,m}^e - \frac{n - m - 1}{n + m}*(r'/R_ref)^2*c_{n-2,m}^e
//
//      -------------------------------------------------------------------
//
//      (1) Basis function          : \bar{c}_{0,0}^e   = [1; 0]
//
//      (2) Diagonal recurrences    : \bar{c}_{n,n}^e   = (2n - 1)*\sqrt{ (1 + \delta_{1,n})/( (2n)*(2n + 1) ) }*[x'/R_ref, -y'/R_ref; y'/R_ref, x'/R_ref]*\bar{c}_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : \bar{c}_{n,n-1}^e = \frac{(2n - 1)}{\sqrt{2n + 1}}*(z'/R_ref)*\bar{c}_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : \bar{c}_{n,m}^e   = (2n - 1)*\sqrt{\frac{2n - 1}{(2n + 1)*(n^2 - m^2)}}*(z'/R_ref)*\bar{c}_{n-1,m}^e - \sqrt{\frac{(2n - 3)*( (n - 1)^2 - m^2 )}{(2n + 1)*(n^2 - m^2)}}*(r'/R_ref)^2*\bar{c}_{n-2,m}^e
//
//      -------------------------------------------------------------------
//
//      where subdiagonal recurrences are determined by m = n - 1 through the vertical recurrences.
//
//      -------------------------------------------------------------------
//
// Note:
//
//  - Modified from the original version 
//
// Dependencies:
//
//  - None
//
// Call
//
//  - None
//
// Called by
//
// - TBD
//
// Modification History:
//  3/22/2015   Siamak Hesar    Modified from the original code written by Yu Takahashi (AccelInteriorPotential_mex.c)
//                              to compute the normalized accelerations and full partials matrix for an exterior gravity field.
// 
//  3/22/2015   Siamak Hesar    Added statements for validating the mex function input types and sizes.
//  06/09/2019  Vishal Ray      Modified function to take order as an input for non-square fields and removed partials w.r.t C and S to increase speed
//////////////////////////////////////////////////////////////////////////*/

#include    "mex.h"
#include	<math.h>
#include	<stdarg.h>
#include	<stdio.h>
#include	<stdlib.h>
#include    <string.h>

#define ABS(x) ((x) < 0) ? -(x) : (x)

#define G 6.67384E-20

#define n_max  2190

#define num_C_max  1325
#define num_S_max  1275

/*//////////////////
// -- Outputs -- //
//////////////////*/

double *Accel_ptr, *Jacobian_ptr;

/*/////////////////
// -- Inputs -- //
/////////////////*/

int    n_degree, m_order;

double R_ref, mu;
double *r_vec_ptr, x_sat, y_sat, z_sat, r_sat;
double *Cbar_ptr, *Sbar_ptr;
double Cbar[n_max + 1][n_max + 1], Sbar[n_max + 1][n_max + 1];

/*////////////////
// -- Index -- //
////////////////*/

int ii, jj, mm, nn, kk;
int Index_C, Index_S;  

/*///////////////////////////
// -- Degree and Order -- //
///////////////////////////*/

double n, m;
double delta_0_m, delta_1_n, delta_1_m, delta_2_m;
double g_bar_int_1, g_bar_int_2, g_bar_int_3;
double s_bar_int_1, s_bar_int_2, s_bar_int_3, s_bar_int_4, s_bar_int_5, s_bar_int_6;

/*///////////////////////////////
// -- Satellite Parameters -- //
///////////////////////////////*/

double x_ddot, y_ddot, z_ddot;
double ddU_dxdx, ddU_dxdy, ddU_dxdz;
double ddU_dydy, ddU_dydz, ddU_dzdz;
double K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10;

double dxddot_dx_bbar_1, dxddot_dy_bbar_1, dxddot_dz_bbar_1, dyddot_dz_bbar_1, dzddot_dz_bbar_1;
double dxddot_dC_bbar_1[num_C_max][num_C_max], dxddot_dS_bbar_1[num_S_max][num_S_max];
double dyddot_dC_bbar_1[num_C_max][num_C_max], dyddot_dS_bbar_1[num_S_max][num_S_max];
double dzddot_dC_bbar_1[num_C_max][num_C_max], dzddot_dS_bbar_1[num_S_max][num_S_max];

/*/////////////////
// -- Output -- //
/////////////////*/

double b_bar_real[n_max + 3][n_max + 3] = {0.0};
double b_bar_imag[n_max + 3][n_max + 3] = {0.0};

/*////////////////////
// -- Functions -- //
////////////////////*/

void GetBnmNormalizedExterior(void);

/***************************************************************
 * main program
 ***************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    /*////////////////////////////////////////
    // Added by Siamak Hesar ////////////////
    // Validating the inputs and outputs ////
    ////////////////////////////////////////*/
    
    /* make sure there are six inputs to the function */
    if(nrhs!=7) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "Seven inputs required.");
    }
    
    /* make sure there are four outputs to the function */
    if(nlhs!=2) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "Two output required.");
    }
      
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                          "First input must be a scalar.");
    }
    
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                          "Second input must be a scalar.");
    }    
    /* make sure the third input argument is scalar */
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                          "Third input must be a scalar.");
    }

    /* make sure the fourth input argument is scalar */
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3]) ||
         mxGetNumberOfElements(prhs[3])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                          "Fourth input must be a scalar.");
    }

    /* make sure the fifth input argument is a row vector of 3 components of type double*/
    if( !mxIsDouble(prhs[4]) || 
         mxIsComplex(prhs[4]) ||
         mxGetM(prhs[4])!=1 ||
         mxGetN(prhs[4])!=3) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Fifth input must be a row vector of 3 components of type double.");
    }
    
    /* make sure the sixth input argument is of type double*/
    if( !mxIsDouble(prhs[5]) || 
         mxIsComplex(prhs[5])) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Sixth input matrix must be of type double.");
    }

    /* make sure the seventh input argument is of type double*/
    if( !mxIsDouble(prhs[6]) || 
         mxIsComplex(prhs[6])) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Seventh input matrix must be of type double.");
    }
   
    /*/////////////////
    // -- Inputs -- //
    /////////////////*/

    n_degree  = (int) mxGetScalar(prhs[0]); /* [n.d.] Degree   */
    m_order   = (int) mxGetScalar(prhs[1]); /* [n.d.] Order   */
    R_ref     = mxGetScalar(prhs[2]); 		/* [km]   Reference Radius*/
    mu	      = mxGetScalar(prhs[3]); 		/* [kg]   Reference mass*/
    r_vec_ptr = mxGetPr(prhs[4]);     		/* [km]   Field point position vector*/
    Cbar_ptr  = mxGetPr(prhs[5]);     		/* [n.d.] C normalized sherical harmonics*/
    Sbar_ptr  = mxGetPr(prhs[6]);     		/* [n.d.] S normalized sherical harmonics*/
    
    /* Added By Siamak Hesar
     * Make sure the C and S matrix sizes are correct*/
     
    if(  mxGetM(prhs[5])!= n_degree + 1 ||
         mxGetN(prhs[5])!= n_degree + 1) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector",
                           "Sixth input matrix must be of size (n_degree + 1)x(n_degree + 1).");
    }
    
    if(  mxGetM(prhs[6])!= n_degree + 1 ||
         mxGetN(prhs[6])!= n_degree + 1) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector",
                           "Seventh input matrix must be of size (n_degree + 1)x(n_degree + 1).");
    }    
    /* END OF INPUT SIZE CHECK */
          

    
    /*//////////////////////
    // -- Allocation -- //
    //////////////////////*/
    
    /*/ - Satellite position*/
    
    x_sat = r_vec_ptr[0]; y_sat = r_vec_ptr[1]; z_sat = r_vec_ptr[2]; /*// [km] x, y, z position vector of the spacecraft */
    r_sat = sqrt(x_sat*x_sat + y_sat*y_sat + z_sat*z_sat);            /*// [km] Norm of the position vector */
    
    /*// - Spherical harmonics*/
        
    for (mm = 0; mm <= n_degree; mm++) {
        
        for (nn = 0; nn <= n_degree; nn++) {
                     
            Cbar[nn][mm] = Cbar_ptr[(n_degree+1)*mm+nn];
            Sbar[nn][mm] = Sbar_ptr[(n_degree+1)*mm+nn];
            
        } /*// for mm*/
        
    } /*// for nn*/
        
    /*///////////////////
    // -- Outputs -- //
    //////////////////*/
    
    plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);     /*// [km/s^2]   Acceleration by \bar{b}_{n,m}^e basis \bar{K}_{nm} */
    plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);     /*// [1/s^2]    STM by \bar{b}_{nm}^e */


    /*// -- Acceleration*/
    
    Accel_ptr   = mxGetPr(plhs[0]);
    
    /*// -- STM*/
    
    Jacobian_ptr    = mxGetPr(plhs[1]);
    
    /*// -- C STM*/
    
   // Partial_C_ptr  = mxGetPr(plhs[2]);

    /*// -- S STM*/
    
    //Partial_S_ptr  = mxGetPr(plhs[3]);
       
    for (nn = 0; nn < 3; nn++) {
        
        Accel_ptr[nn] = 0.0;
        
        for (mm = 0; mm < 3; mm++) {
            
            Jacobian_ptr[3*mm + nn] = 0.0;
            
        } /*// For mm*/
        
    } /*// For nn*/
    

        
    
    GetBnmNormalizedExterior();
           
    /*/////////////////////////
    // -- Pre-allocation -- //
    /////////////////////////*/
    
    x_ddot = 0.0; y_ddot = 0.0; z_ddot = 0.0; 
    
	ddU_dxdx = 0.0;     ddU_dxdy = 0.0;     ddU_dxdz = 0.0;
	ddU_dydy = 0.0;     ddU_dydz = 0.0;     ddU_dzdz = 0.0;
    
	/*/////////////////////////////////////
	// First Partials of the Potential //
	/////////////////////////////////////*/
	
	K0 = 0.5 * mu / R_ref / R_ref;
	
	for (nn = 0; nn<=n_degree; nn++){
    
		n = (double) nn;
		
        if (nn <= m_order) {
         kk = nn;   
        }
        else{
            kk = m_order;
        }
        
		for (mm = 0; mm<=kk; mm++){
			
			m = (double) mm;
			
			if (mm == 1){
				delta_1_m = 1.0;
			}
			else{
				delta_1_m = 0.0;
			} /*// End of if mm == 1*/
			
			K1 = sqrt( (n+2.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+3.0) );
			K2 = sqrt( (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) );
			K3 = sqrt( 2.0 * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_1_m) / (2.0*n+3.0) );
			
			if (mm == 0){
				
				x_ddot -= 2.0*K0 * ( Cbar[nn][mm]*K1*b_bar_real[nn+1][mm+1] );
				y_ddot -= 2.0*K0 * ( Cbar[nn][mm]*K1*b_bar_imag[nn+1][mm+1] );
				z_ddot -= 2.0*K0 * ( Cbar[nn][mm]*sqrt((n-m+1.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0))*b_bar_real[nn+1][mm] );
				
			}
			else{
				
				x_ddot += K0 * ( -Cbar[nn][mm]*K2*b_bar_real[nn+1][mm+1] -Sbar[nn][mm]*K2*b_bar_imag[nn+1][mm+1] +Cbar[nn][mm]*K3*b_bar_real[nn+1][mm-1] +Sbar[nn][mm]*K3*b_bar_imag[nn+1][mm-1]);
				y_ddot += K0 * ( -Cbar[nn][mm]*K2*b_bar_imag[nn+1][mm+1] +Sbar[nn][mm]*K2*b_bar_real[nn+1][mm+1] -Cbar[nn][mm]*K3*b_bar_imag[nn+1][mm-1] +Sbar[nn][mm]*K3*b_bar_real[nn+1][mm-1]);
				z_ddot -= 2.0*K0 * ( Cbar[nn][mm]*sqrt((n-m+1.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0))*b_bar_real[nn+1][mm] +Sbar[nn][mm]*sqrt((n-m+1.0)*(n+m+1.0)*(2*n+1.0)/(2.0*n+3.0))*b_bar_imag[nn+1][mm] );
				
			} /*// End of if mm == 0*/
			
		} /*// End of for mm*/
		
	} /*// End of for nn*/
	
	/*//////////////////////////////////////
	// Second Partials of the Potential //
	//////////////////////////////////////*/
	
	K0 = 0.25 * mu / (R_ref * R_ref * R_ref);

	for (nn = 0; nn<=n_degree; nn++){
		
		n = (double) nn;
		
        if (nn <= m_order) {
         kk = nn;   
        }
        else{
            kk = m_order;
        }
        
		for (mm = 0; mm<=kk; mm++){
			
			m = (double) mm;
			
			if (mm == 1){
				delta_1_m = 1.0;
			}
			else{
				delta_1_m = 0.0;
			} /*// End of if m == 1 */
			
			if (mm == 2){
				delta_2_m = 1.0;
			}
			else{
				delta_2_m = 0.0;
			} /*// End of if m == 2      */
			
			/*// - Common multipliers */
			K1 = sqrt( (n+m+4.0) * (n+m+3.0) * (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+5.0) );
			K2 = sqrt( (n-m+2.0) * (n-m+1.0) * (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+5.0) );
			K3 = sqrt( 2.0 * (n-m+4.0) * (n-m+3.0) * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_2_m) / (2.0*n+5.0) );
			K4 = sqrt( (n+5.0) * (n+4.0) * (n+3.0) * (n+2.0) * (2.0*n+1.0) / (2.0*n+5.0) );
			K5 = sqrt( (n+3.0) * (n+2.0) * (n+1.0) * n * (2.0*n+1.0) / (2.0*n+5.0) );
			K6 = sqrt( (n+4.0) * (n+3.0) * (n+2.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+5.0) );
			K7 = sqrt( (2.0*n+1.0) / (2.0*n+5.0) );
			K8 = sqrt( (n-m+1.0) * (n+m+3.0) * (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+5.0) );
			K9 = sqrt( 2.0 * (n+m+1.0) * (n-m+3.0) * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_1_m) / (2.0*n+5.0) );
			K10= sqrt( (n+3.0) * (n+2.0) * (n+1.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+5.0) );
			
			/*// Partial expressions */
			if (mm == 0){
				
				ddU_dxdx += 2.0 * K0 * ( Cbar[nn][mm]*K6*b_bar_real[nn+2][mm+2] -(n+2)*(n+1)*K7*Cbar[nn][mm]*b_bar_real[nn+2][mm] );
				ddU_dydy -= 2.0 * K0 * ( Cbar[nn][mm]*K6*b_bar_real[nn+2][mm+2] +(n+2)*(n+1)*K7*Cbar[nn][mm]*b_bar_real[nn+2][mm] );
				ddU_dxdy += 2.0 * K0 * ( Cbar[nn][mm]*K6*b_bar_imag[nn+2][mm+2] );
				ddU_dxdz += 4.0 * K0 * ( Cbar[nn][mm]*K10*b_bar_real[nn+2][mm+1] );
				ddU_dydz += 4.0 * K0 * ( Cbar[nn][mm]*K10*b_bar_imag[nn+2][mm+1] );
				ddU_dzdz += 4.0 * K0 * ( Cbar[nn][mm]*K2 *b_bar_real[nn+2][mm] );
				
			}
			else if (mm == 1){
				
				ddU_dxdx += K0 * ( Cbar[nn][mm]*K4*b_bar_real[nn+2][mm+2] +Sbar[nn][mm]*K4*b_bar_imag[nn+2][mm+2] -3.0*K5*Cbar[nn][mm]*b_bar_real[nn+2][mm] -  K5*Sbar[nn][mm]*b_bar_imag[nn+2][mm]);
				ddU_dydy -= K0 * ( Cbar[nn][mm]*K4*b_bar_real[nn+2][mm+2] +Sbar[nn][mm]*K4*b_bar_imag[nn+2][mm+2] +  K5*Cbar[nn][mm]*b_bar_real[nn+2][mm] +3.0*K5*Sbar[nn][mm]*b_bar_imag[nn+2][mm]);
				ddU_dxdy -= K0 * ( Sbar[nn][mm]*K4*b_bar_real[nn+2][mm+2] -Cbar[nn][mm]*K4*b_bar_imag[nn+2][mm+2] +  K5*Sbar[nn][mm]*b_bar_real[nn+2][mm] +  K5*Cbar[nn][mm]*b_bar_imag[nn+2][mm]);
				ddU_dxdz += 2.0 * K0 * ( Cbar[nn][mm]*K8*b_bar_real[nn+2][mm+1] +Sbar[nn][mm]*K8*b_bar_imag[nn+2][mm+1] -  K9*Cbar[nn][mm]*b_bar_real[nn+2][mm-1] -K9*Sbar[nn][mm]*b_bar_imag[nn+2][mm-1] );
				ddU_dydz -= 2.0 * K0 * ( Sbar[nn][mm]*K8*b_bar_real[nn+2][mm+1] -Cbar[nn][mm]*K8*b_bar_imag[nn+2][mm+1] +  K9*Sbar[nn][mm]*b_bar_real[nn+2][mm-1] -K9*Cbar[nn][mm]*b_bar_imag[nn+2][mm-1] );
				ddU_dzdz += 4.0 * K0 * ( Cbar[nn][mm]*K2*b_bar_real[nn+2][mm]   +Sbar[nn][mm]*K2*b_bar_imag[nn+2][mm]);
			}
			else{
				
				ddU_dxdx += K0 * (  Cbar[nn][mm]*K1*b_bar_real[nn+2][mm+2] +Sbar[nn][mm]*K1*b_bar_imag[nn+2][mm+2] -2.0*K2*Cbar[nn][mm]*b_bar_real[nn+2][mm] -2.0*K2*Sbar[nn][mm]*b_bar_imag[nn+2][mm] +K3*Cbar[nn][mm]*b_bar_real[nn+2][mm-2] +K3*Sbar[nn][mm]*b_bar_imag[nn+2][mm-2]);
				ddU_dydy -= K0 * (  Cbar[nn][mm]*K1*b_bar_real[nn+2][mm+2] +Sbar[nn][mm]*K1*b_bar_imag[nn+2][mm+2] +2.0*K2*Cbar[nn][mm]*b_bar_real[nn+2][mm] +2.0*K2*Sbar[nn][mm]*b_bar_imag[nn+2][mm] +K3*Cbar[nn][mm]*b_bar_real[nn+2][mm-2] +K3*Sbar[nn][mm]*b_bar_imag[nn+2][mm-2]);
				ddU_dxdy += K0 * ( -Sbar[nn][mm]*K1*b_bar_real[nn+2][mm+2] +Cbar[nn][mm]*K1*b_bar_imag[nn+2][mm+2] +K3*Sbar[nn][mm]*b_bar_real[nn+2][mm-2] -K3*Cbar[nn][mm]*b_bar_imag[nn+2][mm-2]);
				ddU_dxdz += 2.0 * K0 * ( Cbar[nn][mm]*K8*b_bar_real[nn+2][mm+1] +Sbar[nn][mm]*K8*b_bar_real[nn+2][mm+1] -  K9*Cbar[nn][mm]*b_bar_real[nn+2][mm-1] -K9*Sbar[nn][mm]*b_bar_imag[nn+2][mm-1] );
				ddU_dydz -= 2.0 * K0 * ( Sbar[nn][mm]*K8*b_bar_real[nn+2][mm+1] -Cbar[nn][mm]*K8*b_bar_imag[nn+2][mm+1] +  K9*Sbar[nn][mm]*b_bar_real[nn+2][mm-1] -K9*Cbar[nn][mm]*b_bar_imag[nn+2][mm-1] );
				ddU_dzdz += 4.0 * K0 * ( Cbar[nn][mm]*K2*b_bar_real[nn+2][mm]   +Sbar[nn][mm]*K2*b_bar_imag[nn+2][mm]);
				
			}  /*// End of if mm == 0 */
				   
		} /*// End of for mm */
		
	} /*// End of for nn */

	/*///////////////////////////////
	// Partial Accel partial C&S //
	///////////////////////////////*/
	
	

	/*////////////////////////
	// Assign the outputs //
	////////////////////////*/
	
	Accel_ptr[0] = x_ddot;
	Accel_ptr[1] = y_ddot;
	Accel_ptr[2] = z_ddot;

	Jacobian_ptr[0] = ddU_dxdx;
	Jacobian_ptr[1] = ddU_dxdy;
	Jacobian_ptr[2] = ddU_dxdz;
	Jacobian_ptr[3] = ddU_dxdy;
	Jacobian_ptr[4] = ddU_dydy;
	Jacobian_ptr[5] = ddU_dydz;
	Jacobian_ptr[6] = ddU_dxdz;
	Jacobian_ptr[7] = ddU_dydz;
	Jacobian_ptr[8] = ddU_dzdz;
		
	return;
	
} /*// END of main function*/

/*///////////////////////////////////////////////////////////////////////////*/
void GetBnmNormalizedExterior(void)
{
    
    /*//////////////////////////
    // -- Pre-allocation -- //
    //////////////////////////*/
    
/*    for (nn = 0; nn <= n_degree; nn++) {
        
        for (mm = 0; mm <= n_degree; mm++) {
            
            b_bar_real[nn][mm] = 0.0;
            b_bar_imag[nn][mm] = 0.0;
           
        } // For mm
        
    } // For nn*/
    
    /*////////////////////////////////
    // -- Vertical Recurrences -- //
    ////////////////////////////////*/
    
    for (mm = 0; mm <= n_degree+2; mm++) {
        
        m = (double) mm;
        
        for (nn = mm; nn <= n_degree+2; nn++) {
            
            n = (double) nn;
            
            /*// Recursive Formulae*/
            
            if (mm == nn) {
                
                if (mm == 0) {
                    
                    b_bar_real[0][0] = R_ref/r_sat;
                    b_bar_imag[0][0] = 0.0;
                    
                } else {
                    
                    if (nn == 1) {
                        
                        delta_1_n = 1.0;
                        
                    } else {
                        
                        delta_1_n = 0.0;
                        
                    } /*// For if*/
                    
                    b_bar_real[nn][nn] = sqrt( (1.0 + delta_1_n)*(2.0*n + 1.0)/(2.0*n) ) * (R_ref/r_sat) * ( x_sat/r_sat*b_bar_real[nn-1][nn-1] - y_sat/r_sat*b_bar_imag[nn-1][nn-1] );
                    b_bar_imag[nn][nn] = sqrt( (1.0 + delta_1_n)*(2.0*n + 1.0)/(2.0*n) ) * (R_ref/r_sat) * ( y_sat/r_sat*b_bar_real[nn-1][nn-1] + x_sat/r_sat*b_bar_imag[nn-1][nn-1] );

                } /*// For if*/
                
            } /*// For the Diagonals*/
            
            else {
                
                if ( nn >= 2 ) {
                    
                    b_bar_real[nn][mm] = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_real[nn-1][mm] - sqrt( (2.0*n + 1.0)*( (n - 1.0)*(n - 1.0) - m*m )/( (2.0*n - 3.0)*(n*n - m*m) ) )*(R_ref/r_sat)*(R_ref/r_sat)*b_bar_real[nn-2][mm];
                    b_bar_imag[nn][mm] = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_imag[nn-1][mm] - sqrt( (2.0*n + 1.0)*( (n - 1.0)*(n - 1.0) - m*m )/( (2.0*n - 3.0)*(n*n - m*m) ) )*(R_ref/r_sat)*(R_ref/r_sat)*b_bar_imag[nn-2][mm];
                     
                } else {
                    
                    b_bar_real[nn][mm] = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_real[nn-1][mm];
                    b_bar_imag[nn][mm] = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_imag[nn-1][mm];
                    
                } /*// For if */
                
            } /*// For the Verticals */
            
        } /*// For nn */
        
    } /*// For mm */
        
    return;
    
} /*// END of GetBnmNormalizedExterior */

/*///////////////////////////////////////////////////////////////////////////*/