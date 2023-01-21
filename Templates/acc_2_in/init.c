#include <float.h>
#include "pluto.h"

#define NUMBER_OF_ANALYSIS_VALUES 8
#define AN_MASS                   0
#define AN_E                      1
#define AN_PERI                   2
#define AN_MASS_INNER             3
#define AN_E_INNER                4
#define AN_PERI_INNER             5
#define AN_SIGMA_MIN              6
#define AN_ACC                    7

////////////////////////////////////
// Make these global variables
double g_dm_planet1 = 0.0;
double g_dm_planet2 = 0.0;
////////////////////////////////////

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
    static int first_call = 1;
    g_isoSoundSpeed = g_inputParam[ASPECT_RATIO];
    g_smallDensity = g_inputParam[SIGMA_FLOOR]
                    *g_inputParam[SIGMA_REF]
                    *pow(g_inputParam[A_BIN], -g_inputParam[ALPHA_SIGMA]);

    if (first_call)
    {
        int input_var[256];

        for (int k = 0; k < 256; k++)
            input_var[k] = -1;

        input_var[0] = RHO;
        input_var[1] = VX1;
        input_var[2] = VX2;

        InputDataSet("./in/grid.out", input_var);
        InputDataRead("./in/disc.dbl", " ");
        first_call = 0;
    }
    InputDataInterpolate(v,x1,x2,x3);
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
    double analysisValues[NUMBER_OF_ANALYSIS_VALUES] = {0.0};
    analysisValues[AN_SIGMA_MIN] = DBL_MAX;
    Operators reductionOperators[NUMBER_OF_ANALYSIS_VALUES] = {
        OP_SUM, 
        OP_SUM, 
        OP_SUM, 
        OP_SUM, 
        OP_SUM, 
        OP_SUM, 
        OP_MIN,
        OP_SUM, // AN_ACC
    };

    #if USE_CUDA_ANALYSIS == YES
    callAnalysisReductionKernel(analysisValues, reductionOperators);
    #else
    int i, j, k;
    KDOM_LOOP(k)
    {
        JDOM_LOOP(j)
        {
            double sin_phi = sin(grid[JDIR].x[j]);
            double cos_phi = cos(grid[JDIR].x[j]);

            IDOM_LOOP(i)
            {
                double R = grid[IDIR].x[i];
                double dV = grid[IDIR].dV[i] 
                           *grid[JDIR].dV[j]
                           *grid[KDIR].dV[k];

                double sigma = d->Vc[RHO][k][j][i];
                double v_R = d->Vc[VX1][k][j][i];
                double v_phi = d->Vc[VX2][k][j][i];

                double dm = sigma*dV;

                double factor_1 = -R*v_R*v_phi
                                  /(CONST_G_CODE_UNITS*g_inputParam[M_CO]);
                double factor_2 =  R*v_phi*v_phi
                                  /(CONST_G_CODE_UNITS*g_inputParam[M_CO]) - 1.0;

                double e_x = -factor_1*sin_phi + factor_2*cos_phi;
                double e_y =  factor_1*cos_phi + factor_2*sin_phi; 

                double e_cell = sqrt(e_x*e_x + e_y*e_y);

                double peri_cell = atan2(e_y, e_x);
                if (peri_cell < 0.0)
                    peri_cell += 2*CONST_PI;

                analysisValues[AN_MASS] += dm;
                analysisValues[AN_E] += dm*e_cell;
                analysisValues[AN_PERI] += dm*peri_cell;

                if (R <= g_inputParam[CUTOFF_RADIUS])
                {
                    analysisValues[AN_MASS_INNER] += dm;
                    analysisValues[AN_E_INNER] += dm*e_cell;
                    analysisValues[AN_PERI_INNER] += dm*peri_cell;
                }
            
                analysisValues[AN_SIGMA_MIN] = MIN(analysisValues[AN_SIGMA_MIN],
                                                   sigma);
            }
        }
    }
    #endif

    #ifdef PARALLEL
    double tmp;
    for (int n = 0; n < NUMBER_OF_ANALYSIS_VALUES; n++)
    {
        switch(reductionOperators[n]) {
            case OP_SUM: MPI_Allreduce(analysisValues+n, 
                                       &tmp, 
                                       1, 
                                       MPI_DOUBLE, 
                                       MPI_SUM, 
                                       MPI_COMM_WORLD);
                         analysisValues[n] = tmp;
                         break;
            case OP_PROD: MPI_Allreduce(analysisValues+n, 
                                        &tmp, 
                                        1, 
                                        MPI_DOUBLE, 
                                        MPI_PROD, 
                                        MPI_COMM_WORLD);
                         analysisValues[n] = tmp;
                         break;
            case OP_MIN: MPI_Allreduce(analysisValues+n, 
                                       &tmp, 
                                       1, 
                                       MPI_DOUBLE, 
                                       MPI_MIN, 
                                       MPI_COMM_WORLD);
                         analysisValues[n] = tmp;
                         break;
            case OP_MAX: MPI_Allreduce(analysisValues+n, 
                                       &tmp, 
                                       1, 
                                       MPI_DOUBLE, 
                                       MPI_MAX, 
                                       MPI_COMM_WORLD);
                         analysisValues[n] = tmp;
                         break;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    analysisValues[AN_E]  /= analysisValues[AN_MASS];
    analysisValues[AN_PERI]  /= analysisValues[AN_MASS];

    analysisValues[AN_E_INNER]  /= analysisValues[AN_MASS_INNER];
    analysisValues[AN_PERI_INNER]  /= analysisValues[AN_MASS_INNER];

    ////////////////////////////////////
    g_dm_planet1 = analysisValues[AN_ACC]*g_dt; // this needs to be a mass, or mass rate
    ////////////////////////////////////
    
    /* Write ascii file "averages.dat" to disk */
    if (prank == 0)
    {
        char fname[512];
        static double tpos = -1.0;
        FILE *fp;

        sprintf(fname, "%s/averages.dat", RuntimeGet()->output_dir);

        /* Write header */
        if (g_stepNumber == 0)
        {
            fp = fopen(fname, "w");
            fprintf(fp, "#  t: Time in code units\n");
            fprintf(fp, "#  m: Total disk mass in code units\n");
            fprintf(fp, "#  e_d: Disk eccentricity\n");
            fprintf(fp, "#  peri: Argument of periapsis\n");
            fprintf(fp, "#  e_d_inner: Disk eccentricity of inner disk\n");
            fprintf(fp, "#  peri_inner: Argument of periapsis of inner disk\n");
            fprintf(fp, "#  sigma_min: Minimum surface density in code units\n");
            fprintf(fp, "#  acc: Accreted [???]\n");
            fprintf(fp, "#%-11s  %-12s  %-12s  %-12s  %-12s  %-12s  %-12s  %-12s %-12s\n", 
                    "t", "m", "e_d", "peri", "e_d_inner", "peri_inner", 
                    "sigma_min", "iter", "acc");
        }
        else 
        {
            if (tpos < 0.0)
            {
                char sline[512];
                fp = fopen(fname, "r");
                while (fgets(sline, 512, fp))
                    ;
                sscanf(sline, "%lf\n", &tpos);
                fclose(fp);
            }
            fp = fopen(fname, "a");
        }

        if (g_time > tpos)
        {
            fprintf(fp, "%-12.6e  %-12.6e  %-12.6e  %-12.6e  %-12.6e  %-12.6e  %-12.6e  %-12ld  %-12.6e\n", 
                    g_time, 
                    analysisValues[AN_MASS],
                    analysisValues[AN_E],
                    analysisValues[AN_PERI],
                    analysisValues[AN_E_INNER],
                    analysisValues[AN_PERI_INNER],
                    analysisValues[AN_SIGMA_MIN],
                    g_stepNumber,
                    analysisValues[AN_ACC]);
        }
        fclose(fp);
    }
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
    int i, j, k;

    if (side == 0)
    {
        /* Density floor */
        DOM_LOOP(k,j,i) 
        {
            if (d->Vc[RHO][k][j][i] < g_smallDensity)
            {
                d->Vc[RHO][k][j][i] = g_smallDensity;
            }
        }

        /* Wave damping */
        double Rmax = g_domEnd[IDIR];
        /* damping_factor * (orbital period at R_max) */
        double tau = g_inputParam[DAMPING_FACTOR]
                    *2*CONST_PI
                    *sqrt(Rmax*Rmax*Rmax/(CONST_G_CODE_UNITS*g_inputParam[M_CO]));

        DOM_LOOP(k,j,i)
        {
            double R = grid[IDIR].x[i];
            double sigma_0 = g_inputParam[SIGMA_REF] 
                            *pow(R, -g_inputParam[ALPHA_SIGMA]);

            if (R >= 0.9*Rmax)
            {
                double RoverRmax = R/Rmax;
                double ramp = 100.0*RoverRmax*RoverRmax - 180*RoverRmax + 81.0;
                double lambda = g_dt/tau * ramp;

                d->Vc[RHO][k][j][i] -= lambda * (d->Vc[RHO][k][j][i] - sigma_0);
                /* vR_0 = 0.0 */
                d->Vc[VX1][k][j][i] -= lambda * (d->Vc[VX1][k][j][i]);
            }
        }
    } 
    else if (side == X1_BEG)
    {
        X1_BEG_LOOP(k, j, i)
        {
            /* drho /dr = 0 */
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][2*IBEG-i-1];

            if (d->Vc[VX1][k][j][2*IBEG-i-1] > 0.0)
                d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG-i-1];
            else
                d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][j][2*IBEG-i-1];

            /* domega / dr = 0 */
            d->Vc[VX2][k][j][i] = grid[IDIR].x[i] / grid[IDIR].x[2*IBEG-i-1]
                                 *d->Vc[VX2][k][j][2*IBEG-i-1];
        }
    }
    else if (side == X1_END)
    {
        X1_END_LOOP(k, j, i)
        {
            double R = grid[IDIR].x[i];
            d->Vc[RHO][k][j][i] =  g_inputParam[SIGMA_REF]
                                  *pow(R, -g_inputParam[ALPHA_SIGMA]);
            d->Vc[VX1][k][j][i] = 0.0;
            d->Vc[VX2][k][j][i] = sqrt(CONST_G_CODE_UNITS*g_inputParam[M_CO]/R);
        }
    }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
    g[IDIR] = 0.0;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
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
