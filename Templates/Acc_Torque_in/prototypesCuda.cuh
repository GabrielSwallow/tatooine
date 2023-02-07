/* --------------------------------------------------------------
    Definiton of __device__ function prototypes
   -------------------------------------------------------------- */
#ifndef PROTOTYPES_CUDA_CUH
#define PROTOTYPES_CUDA_CUH

__device__ real hllSolver(int dir, int stateIdx, int stateIdxP1, 
                          int fluxIdx, uint3 dataIdx);

__device__ unsigned int checkZone(uint3 idx, int dir, int offset, 
                                  unsigned int flag);

/* Prototypes from initCuda.cu */
__device__ void cudaUserDefBoundary(int side, uint3 dataIdx, dim3 dataDim);

#if BODY_FORCE != NO
__device__ real cudaBodyForceVector(int dir, uint3 dataIdx, dim3 dataDim);

#if ADD_COOLING_SOURCE == YES
//__device__ real cudaBodyForcePotential(real *v, real x1, real x2, real x3);
__device__ real cudaBodyForcePotential(real x1, real x2, real x3);
#else
__device__ real cudaBodyForcePotential(real x1, real x2, real x3);
#endif
#endif

#if NBODY_SYS == YES
__device__ real cudaNbodySmoothingSquared(int l, uint3 dataIdx, dim3 dataDim);
__device__ real cudaNbodyCalcCellAcceleration(int dir, uint3 dataIdx, dim3 dataDim);
__device__ real cudaNbodyAccretion(int l, uint3 dataIdx, dim3 dataDim);
#endif

#if VISCOSITY != NO
__device__ void cudaViscNu(real *v, real x1, real x2, real x3,
				                    real *nu1, real *nu2);
#endif

#if THERMAL_CONDUCTION != NO
__device__ void cudaTcKappa(real *v, real x1, real x2, real x3,
                            real *kpar, real *knor, real *phi);
#endif

#if USE_CUDA_REDUCTION == YES
__device__ void calcReductionValues(real *values, uint3 dataIdx, dim3 dataDim);
#endif

#if USE_CUDA_ANALYSIS == YES
__device__ void calcAnalysisValues(real *values, uint3 dataIdx, dim3 dataDim);
#endif

#if USE_CPU_CALC_PARAMETERS == YES
void calcCpuParams(real *values);
#endif

#if ADD_COOLING_SOURCE == YES
__device__ real cudaCoolingSource(uint3 dataIdx, dim3 dataDim);
#endif

#endif /* PROTOTYPES_CUDA_CUH */
