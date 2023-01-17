extern "C" {

#include "pluto.h"
#include "plutoCuda.cuh"

/* ********************************************************************************* */
__device__ void cudaUserDefBoundary (int side, uint3 dataIdx, dim3 dataDim)
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 * \param [in] dataIdx
 * \param [in] dataDim
 *
 * The primitive variables of the ghost cells can be set the following way:
 *      cudaV.rho[gpu1D(dataIdx, dataDim)] = ... 
 *      cudaV.vx1[gpu1D(dataIdx, dataDim)] = ... 
 *      cudaV.vx2[gpu1D(dataIdx, dataDim)] = ... (only if COMPONENTS > 1)
 *      cudaV.vx3[gpu1D(dataIdx, dataDim)] = ... (only if COMPONENTS > 2)
 *      cudaV.prs[gpu1D(dataIdx, dataDim)] = ... 
 *
 * The coordinates of the cell can be accessed the following way:
 *      x1 coordinate: cudaGrid.x1[dataIdx.x];
 *      x2 coordinate: cudaGrid.x2[dataIdx.y];
 *      x3 coordinate: cudaGrid.x3[dataIdx.z];
 *
 *********************************************************************************** */
{
    int idx = gpu1D(dataIdx, dataDim);
    real R = cudaGrid.x1[dataIdx.x];
    real phi = cudaGrid.x2[dataIdx.y];
    real sigma_0 = cudaInputParam[SIGMA_REF]
                  *pow(R, -cudaInputParam[ALPHA_SIGMA]);
    
    real xgrid = R * cos(phi);
    real ygrid = R * sin(phi);


    if (side == 0)
    {   
        /* Density floor */
        if (cudaV.rho[idx] < cudaSmallDensity)
        {
            cudaV.rho[idx] = cudaSmallDensity;
        }

        ////////////////////////////////////
        /* Accretion by N-bodies */
        for (int l = 0; l < NB_N; l++)
        {
            real delta_x = xgrid - cudaNb.x[l];
            real delta_y = ygrid - cudaNb.y[l];
            real hill_radius = 0.075;

            // TODO: make this more physical
            real accretion_time_scale = 1/0.069;

            if (delta_x*delta_x + delta_y*delta_y <= hill_radius * hill_radius)
            {
                double acc = cudaV.rho[idx] * cuda_dt / accretion_time_scale;
                // g_dm_planet1 += acc;
                cudaV.rho[idx] -= acc;
            }
        }
        ////////////////////////////////////
        
        
        /* Wave damping */
        real Rmax = cudaDomEnd[IDIR];
        real RoverRmax = R/Rmax;

        /* damping_factor * (orbital period at R_max) */
        real tau = cudaInputParam[DAMPING_FACTOR]
                  *2*CONST_PI
                  *sqrt(Rmax*Rmax*Rmax/(CONST_G_CODE_UNITS*cudaInputParam[M_CO]));

        // if (R >= 0.9*Rmax)
        if (R >= 0.9*Rmax)
        {
            real ramp = 100.0*RoverRmax*RoverRmax - 180*RoverRmax + 81.0;
            real lambda = cuda_dt/tau * ramp;

            cudaV.rho[idx] -= lambda * (cudaV.rho[idx] - sigma_0);
            /* vR_0 = 0.0 */
            cudaV.vx1[idx] -= lambda * (cudaV.vx1[idx]);
        }
    }
    else if (side == X1_BEG)
    {
        int i_act = 2*cudaIBEG - dataIdx.x - 1;
        int offsetIdx = dataIdx.z * dataDim.x * dataDim.y
                       +dataIdx.y * dataDim.x
                       +i_act;

        /* drho / dr = 0 */
        cudaV.vc[RHO][idx] = cudaV.vc[RHO][offsetIdx];

        if (cudaV.vc[VX1][offsetIdx] > 0.0)
            cudaV.vc[VX1][idx] = -cudaV.vc[VX1][offsetIdx];
        else
            cudaV.vc[VX1][idx] = cudaV.vc[VX1][offsetIdx];

        /* domega / dr = 0 */
        cudaV.vc[VX2][idx] =  cudaGrid.x1[dataIdx.x] /cudaGrid.x1[i_act]
                             *cudaV.vc[VX2][offsetIdx];
    }
    else if (side == X1_END)
    {
        cudaV.vc[RHO][idx] = sigma_0;
        cudaV.vc[VX1][idx] = 0.0;
        cudaV.vc[VX2][idx] = sqrt(CONST_G_CODE_UNITS *cudaInputParam[M_CO]/R);
    }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
__device__ real cudaBodyForceVector(int dir, uint3 dataIdx, dim3 dataDim)
/*!
 * Returns the component of the acceleration vector in a given direction.
 * 
 * \param [in] dir      specifies the component of the acceleration vector.
 * \param [in] dataIdx 
 * \param [in] dataDim  
 *
 * The coordinates of the cell can be accessed the following way:
 *      x1 coordinate: cudaGrid.x1[dataIdx.x];
 *      x2 coordinate: cudaGrid.x2[dataIdx.y];
 *      x3 coordinate: cudaGrid.x3[dataIdx.z];
 *
 * The cell-centered primitive variables of the cell can be accessed
 * the following way:
 *      rho: cudaV.rho[gpu1D(dataIdx, dataDim)]
 *      vx1: cudaV.vx1[gpu1D(dataIdx, dataDim)]  
 *      vx2: cudaV.vx2[gpu1D(dataIdx, dataDim)] (only if COMPONENTS > 1)
 *      vx3: cudaV.vx3[gpu1D(dataIdx, dataDim)] (only if COMPONENTS > 2)
 *      prs: cudaV.prs[gpu1D(dataIdx, dataDim)]
 *
 *********************************************************************** */
{
    if (dir == IDIR)
    {
        return 0.0;
    }
    else if (dir == JDIR)
    {
        return 0.0;
    }
    else /* dir == KDIR */
    {
        return 0.0; 

    }
}

/* ********************************************************************* */
__device__ real cudaBodyForcePotential(real x1, real x2, real x3)
/*!
 * Returns the gravitational potential as a function of the coordinates.
 * 
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
    return 0.0;
}
#endif

#if VISCOSITY != NO
__device__ void cudaViscNu(real *v, double x1, double x2, double x3,
                                    double *nu1, double *nu2)
{
    #if BINARY_H == YES
    real x = x1*cos(x2);
    real y = x1*sin(x2);

    real K = 0.0;
    for (int l = 0; l < CENTRAL_OBJECT; l++)
    {
        real dx = x - cudaNb.x[l]; 
        real dy = y - cudaNb.y[l]; 

        real r3 = dx*dx + dy*dy;
        r3 *= sqrt(r3);

        K += CONST_G_CODE_UNITS*cudaNb.m[l] / r3;
    }
    
    *nu1 =  cudaInputParam[ALPHA_VISC]
           *cudaInputParam[ASPECT_RATIO]*cudaInputParam[ASPECT_RATIO]
           *x1*x1
           *sqrt(K)
           *v[RHO];
    #else 
    *nu1 =  cudaInputParam[ALPHA_VISC]
           *cudaInputParam[ASPECT_RATIO]*cudaInputParam[ASPECT_RATIO]
           *sqrt(CONST_G_CODE_UNITS*cudaInputParam[M_CO]*x1)
           *v[RHO];
    #endif
    *nu2 = 0.0;
}
#endif

#if USE_CUDA_REDUCTION == YES
__device__ void calcReductionValues(real *values, uint3 dataIdx, dim3 dataDim)
{
}
#endif

#if USE_CUDA_ANALYSIS == YES
__device__ void calcAnalysisValues(real *values, uint3 dataIdx, dim3 dataDim)
{
    int idx = gpu1D(dataIdx, dataDim);

    real R = cudaGrid.x1[dataIdx.x];
    real phi = cudaGrid.x2[dataIdx.y];

    real xgrid = R * cos(phi);
    real ygrid = R * sin(phi);

    real sin_phi = sin(cudaGrid.x2[dataIdx.y]);
    real cos_phi = cos(cudaGrid.x2[dataIdx.y]);

    real dV =  cudaGrid.dVx1[dataIdx.x]
              *cudaGrid.dVx2[dataIdx.y]
              *cudaGrid.dVx3[dataIdx.z];

    real sigma = cudaV.rho[idx];
    real v_R = cudaV.vx1[idx];
    real v_phi = cudaV.vx2[idx];

    real dm = sigma * dV;

    real factor_1 = -R*v_R*v_phi
                    /(CONST_G_CODE_UNITS*cudaInputParam[M_CO]);
    real factor_2 =  R*v_phi*v_phi
                    /(CONST_G_CODE_UNITS*cudaInputParam[M_CO]) - 1.0;

    real e_x = -factor_1*sin_phi + factor_2*cos_phi;
    real e_y =  factor_1*cos_phi + factor_2*sin_phi; 

    real e_cell = sqrt(e_x*e_x + e_y*e_y);
    real peri_cell = atan2(e_y, e_x);

    if (peri_cell < 0.0)
        peri_cell += 2*CONST_PI;

    values[AN_MASS] = dm;
    values[AN_E] = dm*e_cell;
    values[AN_PERI] = dm*peri_cell;

    if (R <= cudaInputParam[CUTOFF_RADIUS])
    {
        values[AN_MASS_INNER] = dm;
        values[AN_E_INNER] = dm*e_cell;
        values[AN_PERI_INNER] = dm*peri_cell;
    }
    else 
    {
        values[AN_MASS_INNER] = 0.0;
        values[AN_E_INNER] = 0.0;
        values[AN_PERI_INNER] = 0.0;
    }

    values[AN_SIGMA_MIN] = sigma;

    ////////////////////////////////////
    // output the change in mass
    // TODO: make this more physical
    values[AN_ACC] = 0.0;
    real accretion_time_scale = 1/0.069;
    for (int l = 2; l < NB_N; l++)
        {
            real delta_x = xgrid - cudaNb.x[l];
            real delta_y = ygrid - cudaNb.y[l];
            real hill_radius = 0.075;

            if (delta_x*delta_x + delta_y*delta_y <= hill_radius * hill_radius)
            {
                double acc = cudaV.rho[idx] * cuda_dt / accretion_time_scale;
                // this will be total accretion
                values[AN_ACC] = acc * dV/cuda_dt;
            }
        }
    ////////////////////////////////////
}
#endif

#if USE_CPU_CALC_PARAMETERS == YES
void calcCpuParams(real *values)
{
}
#endif

} /* extern "C" */
