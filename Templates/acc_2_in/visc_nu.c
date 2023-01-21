/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <stdio.h>
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3, 
             double *nu1, double *nu2)
/*! 
 * \brief Calculate first and second viscosity coefficients as functions of data and coordinates      
 *
 *    \param [in]      v  pointer to data array containing cell-centered quantities
 *    \param [in]      x1 real, coordinate value 
 *    \param [in]      x2 real, coordinate value 
 *    \param [in]      x3 real, coordinate value 
 *    \param [in, out] nu1 pointer to first viscous coefficient
 *    \param [in, out] nu2 pointer to second viscous coefficient
 *    \return This function has no return value.
 * ************************************************************************** */

{
    #if BINARY_H == YES
    double x = x1*cos(x2);
    double y = x1*sin(x2);

    double K = 0.0;
    for (int l = 0; l < CENTRAL_OBJECT; l++)
    {
        double dx = x - g_nb.x[l];
        double dy = y - g_nb.y[l];

        double r3 = dx*dx + dy*dy;
        r3 *= sqrt(r3);

        K += CONST_G_CODE_UNITS*g_nb.m[l] / r3;
    }
    *nu1 = g_inputParam[ALPHA_VISC]
          *g_inputParam[ASPECT_RATIO]*g_inputParam[ASPECT_RATIO]
          *x1*x1
          *sqrt(K)
          *v[RHO];
    *nu2 = 0.0;
    #else
    *nu1 =  g_inputParam[ALPHA_VISC]
           *g_inputParam[ASPECT_RATIO]*g_inputParam[ASPECT_RATIO]
           *sqrt(CONST_G_CODE_UNITS*g_inputParam[M_CO]*x1)
           *v[RHO];
    *nu2 = 0.0;
    #endif
}
