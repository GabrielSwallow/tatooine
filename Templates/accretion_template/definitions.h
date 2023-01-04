#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                POLAR
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     15
#define  NBODY_SYS               YES

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               EXPLICIT
#define  ROTATING_FRAME          NO

/* -- nbody declarations -- */

#define CENTRAL_OBJECT        BINARY
#define CO_FEELS_DISK         NO
#define INDIRECT_TERM         YES
#define NO_OF_PLANETS         1
#define PLANET_FORMAT         ORBIT
#define ACCRETE               YES

/* -- user-defined parameters (labels) -- */

#define  M_CO                    0
#define  Q_BIN                   1
#define  A_BIN                   2
#define  E_BIN                   3
#define  I_BIN                   4
#define  OMEGA_BIN               5
#define  PERI_BIN                6
#define  F_BIN                   7
#define  ALPHA_VISC              8
#define  ASPECT_RATIO            9
#define  SIGMA_REF               10
#define  ALPHA_SIGMA             11
#define  CUTOFF_RADIUS           12
#define  DAMPING_FACTOR          13
#define  SIGMA_FLOOR             14

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH             (0.08*CONST_au)
#define  UNIT_DENSITY            (1.22121*CONST_Msun/UNIT_LENGTH/UNIT_LENGTH)
#define  UNIT_VELOCITY           (sqrt(CONST_G*UNIT_DENSITY*UNIT_LENGTH))
#define  BINARY_H                NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    NO
#define  PRINT_TO_FILE       YES
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       YES
#define  LIMITER             VANLEER_LIM
