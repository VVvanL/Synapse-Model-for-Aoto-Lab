/*
 *
 *
 * run make to compile
 *
 *
 *
 *
 * NOTE:
 * - in the context of indexing, you HAVE TO subtract 1
 * - constants have to take the same values as in synapse_fun.m
 */


/* invalid index --- used to indicate that an index is not set */
#define    INVALID_INDEX           -1       /* used to point to an out-of-index value */

/* kinetic schemes */

#define     kinetics_a1            1
#define     kinetics_a4            2
    
#define     KINETICS_NO            2       /* number of kinetic schemes */

/* AMPAR states */

#define     amparState_F            1
#define     amparState_B            2
#define     amparState_D            3
#define     amparState_O            4

#define     AMPARSTATES_NO          4

/* kinetic schemes constants */
#define     MAX_AMPARSTATES_NO      9       /* number of maximal AMPAR states of all the scheme implemented */

#define     amparState_C0_a1       1
#define     amparState_C1_a1       2
#define     amparState_C2_a1       3
#define     amparState_D3_a1       4
#define     amparState_D4_a1       5
#define     amparState_D5_a1       6
#define     amparState_D6_a1       7
#define     amparState_D7_a1       8
#define     amparState_O_a1        9

#define     factor  1

#define     k1_a1        (   factor*1.5E7    *1/1000.0 )
#define     k2_a1        (   factor*4E6      *1/1000.0 )
#define     k5_a1        (   factor*1.9863E7 *1/1000.0 )

#define     amparState_C0_a4       1
#define     amparState_C1_a4       2
#define     amparState_C2_a4       3
#define     amparState_D3_a4       4
#define     amparState_D4_a4       5
#define     amparState_D5_a4       6
#define     amparState_D6_a4       7
#define     amparState_D7_a4       8
#define     amparState_O_a4        9


#define     k1_a4        (   factor*2.0E7    *1.1/1000.0 )
#define     k2_a4        (   factor*4E6      *1.1/1000.0 )
#define     k5_a4        (   factor*1.9863E7 *1.1/1000.0 )

#define    amparState_InCleft_OnPSD                 1  /* position tags */
#define    amparState_InCleft_OutsidePSD            2
#define    amparState_InES                          3
#define    amparState_First_InsidePSD_Now_Outside   4
#define    amparState_First_OutsidePSD_Now_Inside   5
#define    amparState_Free_OnPSD                    6
#define    amparState_Free_Outside                  7
#define    amparState_Bound_OnPSD                   8
#define    amparState_Bound_Outside                 9 
#define    amparState_Desensitized_OnPSD           10
#define    amparState_Desensitized_Outside         11 
#define    amparState_Open_OnPSD                   12
#define    amparState_Open_Outside                 13

#define    AMPARTAGS_NO                            13  /* number of AMPAR position tags */

/* glu states */
#define    gluState_FreeCleft           1       /* glu unbound and in happily diffusing in cleft */
#define    gluState_FreeES              2       /* glu unbound and in happily diffusing in cleft */
#define    gluState_BoundToAMPAR        3       /* bound to AMPAR */
#define    gluState_BoundToTransp       4       /* bound to transporter */
#define    gluState_Unreleased          5       /* not yet released */
#define    gluState_Absorbed            6       /* absorbed---diffused out of cleft */
#define    gluState_PumpedOut           7       /* absorbed---diffused out of cleft     */
#define    gluState_UnbindingFromTransp 8       /* use if glu has just been unbound from transporter, */
#define    gluState_UnbindingFromAMPAR  9       /* use if glu has just been unbound from AMPAR
                                                 * to indicate that the glu has to be re-inserted
                                                 * into diffusion space correctly */
#define    GLUSTATES_NO                 9       /* number of AMPAR states */

/* transporter states */
#define    transpState_unbound                  1
#define    transpState_bound_nextState_unbound  2
#define    transpState_bound_nextState_pumpedIn 3
#define    transpState_pumpedIn                 4


/* pi */
#define    PI          3.1415926535897931   /* what is pi ? */

/* MATLAB integer types */
#define     uint8   unsigned char
#define     uint16  unsigned short
#define     uint32  unsigned int

#define     int8    char
#define     int16   short
#define     int32   int

/* includes */
#include    <math.h>
#include    <string.h>
#include    "mex.h"

#include    <time.h>



void *getMATLABvalue(mxArray *var_ptr) {
    if (var_ptr == NULL) {
        mexErrMsgTxt("Invalid MATLAB var_ptr.\n");
    }
    return(mxGetData(var_ptr));
}

mxArray *getMATLABfield(const char *structName, const char *fieldName) {
    mxArray *struct_ptr, *field_ptr;
    
    struct_ptr = (mxArray*)mexGetVariablePtr("caller", structName);
    if (struct_ptr == NULL) {
        mexPrintf("MATLAB structure %s: ", structName);
        mexErrMsgTxt("MATLAB structure not found.\n");
    }

    field_ptr = mxGetField(struct_ptr, 0, fieldName);
    return(field_ptr);
}

mxArray *getMATLABvar(const char *varName) {
    mxArray *var_ptr = mexGetVariable("caller", varName);
    if (var_ptr == NULL) {
        mexPrintf("MATLAB variable %s: ", varName);
        mexErrMsgTxt("MATLAB variable not found.\n");
    }
    return(var_ptr);
}

void setMATLABvar_ptr(const char *varName, mxArray *var_ptr) {
    if (var_ptr == NULL) {
        mexErrMsgTxt("mxArray* var_ptr not defined. Could not set MATLAB variable.\n");
    } else {
        mexPutVariable("caller", varName, var_ptr);
    }    
    return;    
}

/********************************/
            double V1, V2, S;
            int phase = 0;
        double gaussrand() {
            double X;
            
            if(phase == 0) {
                do {
                    double U1 = (double)rand() / RAND_MAX;
                    double U2 = (double)rand() / RAND_MAX;
                    
                    V1 = 2 * U1 - 1;
                    V2 = 2 * U2 - 1;
                    S = V1 * V1 + V2 * V2;
                } while(S >= 1 || S == 0);
                
                X = V1 * sqrt(-2 * log(S) / S);
            } else
                X = V2 * sqrt(-2 * log(S) / S);
            
            phase = 1 - phase;
            
            return X;
        }
/***********************************/


/* MATLAB function implementation */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* MATLAB variables */
    
    /* time */
    int     timestep        = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "timestep"));
    int     timeSteps       = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "timeSteps"));
    double  timeStepSize    =      *(double*)getMATLABvalue(getMATLABfield("SIM", "timeStepSize")); /* in ms
    
    /* vesicle releases */
    int     vesicleNo       = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "vesicleNo"));
    
    /* synapse geometry */
    double  R_PSD           =      *(double*)getMATLABvalue(getMATLABfield("SIM", "R_PSD"));
    double  R_Cleft         =      *(double*)getMATLABvalue(getMATLABfield("SIM", "R_Cleft"));
    double  R_ES            =      *(double*)getMATLABvalue(getMATLABfield("SIM", "R_ES"));
    double  R_Reservoir     =      *(double*)getMATLABvalue(getMATLABfield("SIM", "R_Reservoir"));
    double  lowerLid        =      *(double*)getMATLABvalue(getMATLABfield("SIM", "lowerLid"));
    double  upperLid        =      *(double*)getMATLABvalue(getMATLABfield("SIM", "upperLid"));
    double  lowerLidES      =      *(double*)getMATLABvalue(getMATLABfield("SIM", "lowerLidES"));
    double  upperLidES      =      *(double*)getMATLABvalue(getMATLABfield("SIM", "upperLidES"));
    
    int     esBinNo         = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "esBinNo"));
    double  esBinSize       =      *(double*)getMATLABvalue(getMATLABfield("SIM", "esBinSize"));
    
    /* glu */
    int     gluNo           = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "gluNo"));
    double  D_glu           =      *(double*)getMATLABvalue(getMATLABfield("SIM", "D_glu"));
    double  dt_glu          =      *(double*)getMATLABvalue(getMATLABvar  (       "dt_glu"));
    
    /* AMPARs */
    int     amparNo         = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "amparNo"));
    double  R_BindingToAMPAR=      *(double*)getMATLABvalue(getMATLABfield("SIM", "R_BindingToAMPAR"));
    double  D_ampar_PSD     =      *(double*)getMATLABvalue(getMATLABfield("SIM", "D_ampar_PSD"));
    double  D_ampar_outside =      *(double*)getMATLABvalue(getMATLABfield("SIM", "D_ampar_outside"));
    double  dt_ampar_PSD    =      *(double*)getMATLABvalue(getMATLABvar  (       "dt_ampar_PSD"));
    double  dt_ampar_outside=      *(double*)getMATLABvalue(getMATLABvar  (       "dt_ampar_outside"));

    int     P_Reflect_Inside =(int)((*(double*)getMATLABvalue(getMATLABfield("SIM","P_Reflect_Inside")) ) *RAND_MAX);
    int     P_Reflect_Outside=(int)((*(double*)getMATLABvalue(getMATLABfield("SIM","P_Reflect_Outside"))) *RAND_MAX);
    
    double  RK_ErrorTol     =      *(double*)getMATLABvalue(getMATLABfield("SIM", "RK_ErrorTol"));  /* error tolerance */
    double  RK_MaxRef       =      *(double*)getMATLABvalue(getMATLABfield("SIM", "RK_MaxRef"));    /* maximal refinement */
    
    /* transporters */
    double  D_transp        =      *(double*)getMATLABvalue(getMATLABfield("SIM", "D_transp"));
    double  dt_transp       =      *(double*)getMATLABvalue(getMATLABvar("dt_transp"));
    double  transpDensity   =      *(double*)getMATLABvalue(getMATLABfield("SIM", "transpDensity")); /* in nm^-2 */
    double  transpBindingTime   =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transpBindingTime")); /* in ms */
    double  transpPumpingTime   =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transpPumpingTime")); /* in ms */
    double  R_BindingToTransp   =  *(double*)getMATLABvalue(getMATLABfield("SIM", "R_BindingToTransp")); /* in nm */
    double  transp2transpDist   =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transp2transpDist")); /* in nm */
    int     transpHor       = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "transpHor"));
    int     transpVer       = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "transpVer"));
    int     transpNo        = (int)*(double*)getMATLABvalue(getMATLABfield("SIM", "transpNo"));
    int     P_bindingToTransp=(int)((*(double*)getMATLABvalue(getMATLABfield("SIM", "P_bindingToTransp"))) *RAND_MAX);

    double  transpRate_assoc    =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transpRate_assoc"));
    double  transpRate_dissoc   =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transpRate_dissoc"));
    double  transpRate_pump     =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transpRate_pump"));
    double  transpRate_getReady =  *(double*)getMATLABvalue(getMATLABfield("SIM", "transpRate_getReady"));
    
    /* kinetics_struct */
    int     *kinetics_amparState_F =   (int*)getMATLABvalue(getMATLABfield("kinetics", "amparState_F")); /*  zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ... */
    int     *kinetics_amparState_B =   (int*)getMATLABvalue(getMATLABfield("kinetics", "amparState_B")); /*  zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ... */
    int     *kinetics_amparState_D =   (int*)getMATLABvalue(getMATLABfield("kinetics", "amparState_D")); /*  zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ... */
    int     *kinetics_amparState_O =   (int*)getMATLABvalue(getMATLABfield("kinetics", "amparState_O")); /*  zeros(KINETICS_NO, MAX_AMPARSTATES_NO, 'int32'), ... */
    double  *kinetics_Q_transition =(double*)getMATLABvalue(getMATLABfield("kinetics", "Q_transition")); /*  zeros(KINETICS_NO, MAX_AMPARSTATES_NO, MAX_AMPARSTATES_NO, 'double'), ... */
    double  *kinetics_P_transition =(double*)getMATLABvalue(getMATLABfield("kinetics", "P_transition")); /*  zeros(KINETICS_NO, MAX_AMPARSTATES_NO, MAX_AMPARSTATES_NO, 'double'), ... */
    
    /* vesicle_struct */
    double  *vesicle_start          =       (double*)getMATLABvalue(getMATLABfield("vesicle", "start"));          /* zeros(1, vesicleNo) */
	double  *vesicle_duration       =       (double*)getMATLABvalue(getMATLABfield("vesicle", "duration"));       /* zeros(1, vesicleNo) */
    double  *vesicle_gluNo          =       (double*)getMATLABvalue(getMATLABfield("vesicle", "gluNo"));          /* zeros(1, vesicleNo) */
    double  *vesicle_gluNo_Released =       (double*)getMATLABvalue(getMATLABfield("vesicle", "gluNo_Released"));          /* zeros(1, vesicleNo) */
    double  *vesicle_position       =       (double*)getMATLABvalue(getMATLABfield("vesicle", "position"));       /* zeros(3, vesicleNo) */
    double  *vesicle_timeStep_Start =       (double*)getMATLABvalue(getMATLABfield("vesicle", "timeStep_Start"));/*zeros(1, vesicleNo) */
    double  *vesicle_timeStep_ReleaseRate = (double*)getMATLABvalue(getMATLABfield("vesicle", "timeStep_ReleaseRate")); /* zeros(1, vesicleNo) */
    double  *vesicle_timeStep_End   =       (double*)getMATLABvalue(getMATLABfield("vesicle", "timeStep_End"));   /* zeros(1, vesicleNo) */

    /* ampar_struct */
    double  *ampar_traj         = (double*)getMATLABvalue(getMATLABfield("ampar", "traj"));             /* zeros(3, 1, amparNo), ... % 3=dimension, 2d diffusion but keep z-coord constant */
    double  *ampar_stateDistr   = (double*)getMATLABvalue(getMATLABfield("ampar", "stateDistr"));       /* zeros(MAX_AMPARSTATES_NO, amparNo, 'double'), ... */
    int8    *ampar_kinetics     =   (int8*)getMATLABvalue(getMATLABfield("ampar", "kinetics"));         /* zeros(1, amparNo, 'int8'), ... */
    uint32  *ampar_currentState = (uint32*)getMATLABvalue(getMATLABfield("ampar", "currentState"));     /* ones (1, amparNo, 'uint32')*INVALID_INDEX, ... */
 /* int32   *ampar_ligand1      =  (int32*)getMATLABvalue(getMATLABfield("ampar", "ligand1"));       */ /* zeros(1, amparNo, 'uint32'), ... */
 /* int32   *ampar_ligand2      =  (int32*)getMATLABvalue(getMATLABfield("ampar", "ligand2"));       */ /* zeros(1, amparNo, 'uint32'), ... */
    int32   *ampar_hits         =  (int32*)getMATLABvalue(getMATLABfield("ampar", "hits"));             /* zeros(1, amparNo, 'uint32'), ... */
    uint8   *ampar_tag          =  (uint8*)getMATLABvalue(getMATLABfield("ampar", "tag"));              /* zeros(1, amparNo, 'uint8'), ... */
    uint8   *ampar_tag_initial  =  (uint8*)getMATLABvalue(getMATLABfield("ampar", "tag_initial"));      /* zeros(1, amparNo, 'uint8'), ... */
    double  *ampar_T            = (double*)getMATLABvalue(getMATLABfield("ampar", "T"));                /* zeros(1, amparNo)); */

    /* transp_struct */
    double  *transp_pos         = (double*)getMATLABvalue(getMATLABfield("transp", "pos"));             /* zeros(3, transpHor, transpVer), ... */
    uint8   *transp_state       =  (uint8*)getMATLABvalue(getMATLABfield("transp", "state"));           /* ones (transpHor, transpVer, 'uint8') .* transpState_unbound, ... */
    uint32  *transp_ligand      = (uint32*)getMATLABvalue(getMATLABfield("transp", "ligand"));          /* zeros(transpHor, transpVer, 'uint32'), ... */
    double  *transp_pumpingTime = (double*)getMATLABvalue(getMATLABfield("transp", "pumpingTime"));     /* zeros(transpHor, transpVer)); */
    double  *transp_T           = (double*)getMATLABvalue(getMATLABfield("transp", "T"));               /* zeros(transpHor, transpVer)); */
    
    /* glu_struct */
    double  *glu_traj           = (double*)getMATLABvalue(getMATLABfield("glu", "traj"));               /* zeros(3, gluNo_AllGlus), ... % 3= dimension, 3d diffusion */
    uint8   *glu_state          =  (uint8*)getMATLABvalue(getMATLABfield("glu", "state"));              /* ones (1, gluNo_AllGlus, 'uint8') .* gluState_Unreleased, ... */
    uint32  *glu_boundToAMPAR   = (uint32*)getMATLABvalue(getMATLABfield("glu", "boundToAMPAR"));       /* zeros(1, gluNo_AllGlus, 'uint32'), ... */
    double  *glu_bindingPosition= (double*)getMATLABvalue(getMATLABfield("glu", "bindingPosition"));    /* zeros(3, gluNo_AllGlus), ... */
    double  *glu_birthTime      = (double*)getMATLABvalue(getMATLABfield("glu", "birthTime"));          /* zeros(1, gluNo_AllGlus), ... */
    double  *glu_exitTime       = (double*)getMATLABvalue(getMATLABfield("glu", "exitTime"));           /* ones (1, gluNo_AllGlus) .* timeSteps); */

    /* general (insertion distribution look-up table) */
    int     insertDistrTableLen = (int)*(double*)getMATLABvalue(getMATLABvar("insertDistrTableLen")) -1; /* subtract 1 because we want to use `rand() AND ...' below */
    double  *insertDistrTable   =   (double*)getMATLABvalue(getMATLABvar("insertDistrTable"));  /* zeros(insertDistrTableLen, 1); */
    
    /* book-keeping variables */
    mxArray *amparTags_ptr, *amparTagsVar_ptr,
            *amparStates_ptr, *amparMean_ptr, *amparVarS_ptr,
            *gluStates_ptr, *gluDistrES_ptr;
    double  *amparTags          = (double*)getMATLABvalue(amparTags_ptr    = getMATLABvar("amparTags"));    /* zeros(timeSteps, AMPARTAGS_NO, 'double'); % array of AMPAR tags */
    double  *amparTagsVar       = (double*)getMATLABvalue(amparTagsVar_ptr = getMATLABvar("amparTagsVar")); /* zeros(timeSteps, AMPARTAGS_NO, 'double'); % array of AMPAR tags */
    double  *amparStates        = (double*)getMATLABvalue(amparStates_ptr  = getMATLABvar("amparStates"));  /* zeros(timeSteps, MAX_AMPARSTATES_NO, 'double'); % array of AMPAR states */
    double  *amparMean          = (double*)getMATLABvalue(amparMean_ptr    = getMATLABvar("amparMean"));    /* zeros(timeSteps, MAX_AMPARSTATES_NO, 'double'); % array of AMPAR states */
    double  *amparVarS          = (double*)getMATLABvalue(amparVarS_ptr    = getMATLABvar("amparVarS"));    /* zeros(timeSteps, MAX_AMPARSTATES_NO, 'double'); % array of AMPAR states */
    uint32  *gluStates          = (uint32*)getMATLABvalue(gluStates_ptr    = getMATLABvar("gluStates"));    /* zeros(timeSteps, GLUSTATES_NO, 'uint32');   % no of glus with same states, 8 states */
    uint32  *gluDistrES         = (uint32*)getMATLABvalue(gluDistrES_ptr   = getMATLABvar("gluDistrES"));   /* zeros(timeSteps, transpVer, 'uint32');     % no of glus in vertical bin-rings in ES in z-direction */
    
    /* local variables */
    int     i_vesicle, i_ampar, i_states, i_glu, i_newly;
    double  X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3;
    double  R, distTimesR;
    double  v[3];
    double  theta;
    int     i_hor, i_ver, checking, checkingES;
    double  hitsToConcentration = sqrt( PI / (D_glu*timeStepSize) ) / (R_BindingToAMPAR*R_BindingToAMPAR*PI) *10.0/6.023;
    
    /* transition probability matrix     */
    double  Q_C1_C1_a1, Q_D3_D3_a1, Q_C1_C1_a4, Q_D3_D3_a4;
    double  Q_transition    [KINETICS_NO][MAX_AMPARSTATES_NO][MAX_AMPARSTATES_NO];
    
    double  P_transition    [KINETICS_NO][MAX_AMPARSTATES_NO][MAX_AMPARSTATES_NO];

    /* RK scheme */
     double  RK_A[6]    =   {          0.0,        1.0/4.0,        3.0/8.0,       12.0/13.0,       1.0,  1.0/2.0};
     
     double  RK_B[6][5] = { {          0.0,            0.0,            0.0,             0.0,       0.0},
                            {      1.0/4.0,            0.0,            0.0,             0.0,       0.0},
                            {     3.0/32.0,       9.0/32.0,            0.0,             0.0,       0.0},
                            {1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0,             0.0,       0.0},
                            {  439.0/216.0,           -8.0,   3680.0/513.0,    845.0/4104.0,       0.0},
                            {    -8.0/27.0,            2.0, -3544.0/2565.0,   1859.0/4104.0, -11.0/40.0} };
        
     double RK_C[2][6]  = { {   16.0/135.0,            0.0, 6656.0/12825.0, 28561.0/56430.0,  -9.0/50.0, 2.0/55.0},
                            {   25.0/216.0,            0.0,  1408.0/2565.0,   2197.0/4104.0,   -1.0/5.0,    0.0} };

    /*>>>> function starts here >>>>>>>>>>>>>>>>>>>> */
    
/*    srand((unsigned)time(NULL));
 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*    srand(0); */
    
    if ((nrhs != 0) || (nlhs != 0)) {
        mexErrMsgTxt("Don't specify input or output arguments.");
        return;
    }

    /* initialize P (transition prob matrix) */
    {
        int i, j, k;
        
        for (k = 0; k < KINETICS_NO; k++) {
            for (i = 0; i < MAX_AMPARSTATES_NO; i++)
                for (j = 0; j < MAX_AMPARSTATES_NO; j++) {
                    Q_transition[k][i][j] = *(kinetics_Q_transition + k + (i + j*MAX_AMPARSTATES_NO)*KINETICS_NO);
                    P_transition[k][i][j] = *(kinetics_P_transition + k + (i + j*MAX_AMPARSTATES_NO)*KINETICS_NO);
                }
        }

    /* initialize diagonal elements that have to be updated due to concentration-dependent rates */
        
    /* case kinetics_a1 : */
            Q_C1_C1_a1 = Q_transition[kinetics_a1-1]   [amparState_C1_a1-1][amparState_C1_a1-1]; /* -1 == MATLAB indexing */
            Q_D3_D3_a1 = Q_transition[kinetics_a1-1]   [amparState_D3_a1-1][amparState_D3_a1-1];
    /* case kinetics_a4: */
            Q_C1_C1_a4 = Q_transition[kinetics_a4-1]   [amparState_C1_a4-1][amparState_C1_a4-1];
            Q_D3_D3_a4 = Q_transition[kinetics_a4-1]   [amparState_D3_a4-1][amparState_D3_a4-1]; /* -1 == MATLAB indexing */
    }
    
/* VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV */
/*  VVVVVVVVV  VVVVVVVVV  VVVVVVVVV  VVVVVVVVV  VVVVVVVVV  VVVVVVVVV  VVVVVVVVV  */
/*   VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV    VVVVVVV   */
/*    VVVVV      VVVVV      VVVVV      VVVVV      VVVVV      VVVVV      VVVVV    */
/*     VVV        VVV        VVV        VVV        VVV        VVV        VVV     */
/*      V          V          V          V          V          V          V      */

/*gluNo = 0; */

timestep = 0; 
while (timestep < timeSteps) {

    /*---------------------------------------------------------------------- */
    /* GLUTAMATE RELEASE */
    /* */
    /*increase timestep AFTER checking how many glu's are to be released */
    for (i_vesicle = 0; i_vesicle < vesicleNo; i_vesicle++) {
        if (timestep >= *(vesicle_timeStep_Start + i_vesicle)  &&
            timestep <= *(vesicle_timeStep_End   + i_vesicle)) {
            double rel = floor( ((double)timestep - *(vesicle_timeStep_Start + i_vesicle)+1) * (*(vesicle_timeStep_ReleaseRate + i_vesicle)) );

            /* can't happen, but still check */
            if (rel > *(vesicle_gluNo + i_vesicle))
                rel = *(vesicle_gluNo + i_vesicle);
            /* in case of numerical problems, release missing glus in last release step */
            if (timestep == *(vesicle_timeStep_End + i_vesicle))
                rel = *(vesicle_gluNo + i_vesicle);            
            
            for (i_newly = *(vesicle_gluNo_Released + i_vesicle); i_newly < rel; i_newly++) {
                *(glu_state + gluNo) = gluState_FreeCleft;
                *(glu_traj + 0 + gluNo*3) = *(vesicle_position + 0 + i_vesicle*3);
                *(glu_traj + 1 + gluNo*3) = *(vesicle_position + 1 + i_vesicle*3);
                *(glu_traj + 2 + gluNo*3) = *(vesicle_position + 2 + i_vesicle*3);
                *(glu_birthTime + gluNo) = timestep;
                gluNo ++;
            }
            *(vesicle_gluNo_Released + i_vesicle) = rel;
                    
        } /* if timestep within release window */
    } /* for i_vesicle */
    
    
    /*---------------------------------------------------------------------- */
    /*---------------------------------------------------------------------- */
    /*---------------------------------------------------------------------- */
    /* WORK HERE */
    /*   1) open channels */
    /*   2) number of free glus */
    /*   3) ... */
    

    /* amparStates */
    /* amparMean --- save 2*amparStates_THISRUN   (where amparStates_THISRUN = ampar_stateDistr) */
    /* amparVarS --- save amparStates_THISRUN ^2 */
    {
        for (i_states = 0; i_states < AMPARSTATES_NO; i_states++) {
            double x_n = 0;            
            int    *amparStates_List_ptr;
            
            switch (i_states) {
                case amparState_F-1 : amparStates_List_ptr = kinetics_amparState_F; break;
                case amparState_B-1 : amparStates_List_ptr = kinetics_amparState_B; break;
                case amparState_D-1 : amparStates_List_ptr = kinetics_amparState_D; break;
                case amparState_O-1 : amparStates_List_ptr = kinetics_amparState_O; break;
            }
            
            for (i_ampar = 0; i_ampar < amparNo; i_ampar++) {
                int k = *(ampar_kinetics + i_ampar) - 1, s;
                
                for (s = 0; s < MAX_AMPARSTATES_NO; s++) {
                    if (*(amparStates_List_ptr + k + s*KINETICS_NO) == 1)
                        x_n += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO);
                }
            /* x_n is the n-th sample of 1..SIM.runs (so n is i_runs) */
            }

            *(amparStates + timestep + i_states *timeSteps) += x_n;
            
            *(amparMean   + timestep + i_states *timeSteps) += x_n;
            *(amparVarS   + timestep + i_states *timeSteps) += x_n * x_n;
     
        }
    }

    /*>>>> gluStates */
    for (i_glu = 0; i_glu < gluNo; i_glu++) {
        (*(gluStates + timestep + ((int)(*(glu_state + i_glu)  -1)) *timeSteps)) ++;
    }
    /*<<<< gluStates */

    /*>>>> amparTags */
    {
        int i_tags, k, s;
        int tags[AMPARTAGS_NO];
  
        memset(tags, 0, AMPARTAGS_NO*sizeof(int));
        
        for (i_ampar = 0; i_ampar < amparNo; i_ampar++) {
            k = *(ampar_kinetics + i_ampar) - 1;
            
            tags[*(ampar_tag + i_ampar) - 1] ++; /* refers only to the first 3 tags */
            
            (*(amparTags + timestep + ((int)(*(ampar_tag + i_ampar)-1)) *timeSteps)) ++;
            
            /* PSD - PSD */
            if (*(ampar_tag_initial + i_ampar) == amparState_InCleft_OnPSD  &&
                *(ampar_tag         + i_ampar) != amparState_InCleft_OnPSD) {
                (*(amparTags + timestep + (amparState_First_InsidePSD_Now_Outside - 1) *timeSteps)) ++;
                tags[amparState_First_InsidePSD_Now_Outside - 1] ++;
            }
            
            if (*(ampar_tag_initial + i_ampar) != amparState_InCleft_OnPSD  &&
                *(ampar_tag         + i_ampar) == amparState_InCleft_OnPSD) {
                (*(amparTags + timestep + (amparState_First_OutsidePSD_Now_Inside - 1) *timeSteps)) ++;
                tags[amparState_First_OutsidePSD_Now_Inside - 1] ++;
            }

            /* report AMPAR tags according to internal state */
            if (*(ampar_tag + i_ampar) == amparState_InCleft_OnPSD) {
                for (s = 0; s < MAX_AMPARSTATES_NO; s++) {
                    *(amparTags + timestep + (amparState_Free_OnPSD         -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_F + k + s*KINETICS_NO));
                    *(amparTags + timestep + (amparState_Bound_OnPSD        -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_B + k + s*KINETICS_NO));
                    *(amparTags + timestep + (amparState_Desensitized_OnPSD -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_D + k + s*KINETICS_NO));
                    *(amparTags + timestep + (amparState_Open_OnPSD         -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_O + k + s*KINETICS_NO));
                }
            } else {
                for (s = 0; s < MAX_AMPARSTATES_NO; s++) {
                    *(amparTags + timestep + (amparState_Free_Outside         -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_F + k + s*KINETICS_NO));
                    *(amparTags + timestep + (amparState_Bound_Outside        -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_B + k + s*KINETICS_NO));
                    *(amparTags + timestep + (amparState_Desensitized_Outside -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_D + k + s*KINETICS_NO));
                    *(amparTags + timestep + (amparState_Open_Outside         -1)*timeSteps) += *(ampar_stateDistr + s + i_ampar *MAX_AMPARSTATES_NO) * (*(kinetics_amparState_O + k + s*KINETICS_NO));
                }
            }
        } /* for i_ampar */
        
        /* amparTagsVar --- variance*/
        for (i_tags = 0; i_tags < AMPARTAGS_NO; i_tags++) {
            double x_n = tags[i_tags];
            (*(amparTagsVar + timestep + i_tags* timeSteps)) += x_n * x_n;
        }
    }
    /*<<<< amparTags */
    
    /*---------------------------------------------------------------------- */
    /*---------------------------------------------------------------------- */
    /*---------------------------------------------------------------------- */
    /* UPDATE: */
    /*   transporters */
    /*   AMPAR states */
    /*   AMPAR diffusion */
    /*   glu diffusion */
    /*   time step */

    /*---------------------------------------------------------------------- */
    /* UPDATE transporters */
    /* */
    /*>>>> transporters */
    if (transpNo > 0) {
        int ligand;
        
        for (i_hor = 0; i_hor < transpHor; i_hor++) {
            for (i_ver = 0; i_ver < transpVer; i_ver++) {
                if (*(transp_state + i_hor + i_ver*transpHor) != transpState_unbound) {

                    /* decrease waiting time */
                    if (*(transp_T + i_hor + i_ver*transpHor) > 0) {
                        *(transp_T + i_hor + i_ver*transpHor) -= timeStepSize;
                        if (*(transp_T + i_hor + i_ver*transpHor) < 0)
                            *(transp_T + i_hor + i_ver*transpHor) = 0;                        
                    } else {
                        /* waiting time ran out */
                        
                        if (*(transp_state + i_hor + i_ver*transpHor) == transpState_bound_nextState_pumpedIn) {
                            /* pump in */
                            *(transp_state + i_hor + i_ver*transpHor) = transpState_pumpedIn;
                            /* deterministic transporters */
                              /* *(transp_T + i_hor + i_ver*transpHor) = transp_getReadyTime; */
                            /* stochastic transporters */
                            *(transp_T + i_hor + i_ver*transpHor) = -log((double)rand()/RAND_MAX) / transpRate_getReady;
                            
                        } else 
                        if (*(transp_state + i_hor + i_ver*transpHor) == transpState_bound_nextState_unbound) {
                            /* unbind and release glu */
                            *(transp_state + i_hor + i_ver*transpHor) = transpState_unbound;
                            ligand = *(transp_ligand + i_hor + i_ver*transpHor);
                            *(transp_ligand + i_hor + i_ver*transpHor) = INVALID_INDEX;

                            /* The correct release position (insertion in the ES) 
                             * will be done in the glu update loop below */
                            *(glu_state + ligand) = gluState_UnbindingFromTransp;
                            /* Set the unbinding position of the glu to the current
                             * position of the transoporter from which it unbinds
                             * ( @@@@ this is for later, when transporter diffusion
                             *  will be included)
                             */
/*
                            *(glu_traj + 0 + ligand*3) = *(transp_pos + 0 + (i_hor + i_ver*transpHor)*3);
                            *(glu_traj + 1 + ligand*3) = *(transp_pos + 1 + (i_hor + i_ver*transpHor)*3);
                            *(glu_traj + 2 + ligand*3) = *(transp_pos + 2 + (i_hor + i_ver*transpHor)*3);
*/
                            X1 = *(transp_pos + 0 + (i_hor + i_ver*transpHor)*3);
                            X2 = *(transp_pos + 1 + (i_hor + i_ver*transpHor)*3);
                            X3 = *(transp_pos + 2 + (i_hor + i_ver*transpHor)*3);

    if (X1*X1 + X2*X2 > R_ES*R_ES) {
        X1 = X1*0.95;
        X2 = X2*0.95;
    }
                            *(glu_traj + 0 + ligand*3) = X1;
                            *(glu_traj + 1 + ligand*3) = X2;
                            *(glu_traj + 2 + ligand*3) = X3;                            
                        } else 
                        if (*(transp_state + i_hor + i_ver*transpHor) == transpState_pumpedIn) {
                            *(transp_state + i_hor + i_ver*transpHor) = transpState_unbound;
                            /* mark glu as pumped-out */
                            ligand = *(transp_ligand + i_hor + i_ver*transpHor);
                            *(transp_ligand + i_hor + i_ver*transpHor) = INVALID_INDEX;
                            *(glu_state + ligand) = gluState_PumpedOut; /* kill glu */
                            *(glu_exitTime + ligand) = timestep*timeStepSize; /* remember time at which pump-out was finished */
                        }
                        
                    }

                } /* if */

            } /* for i_ver */
        } /* for i_hor */
    } /* If transpNo > 0 */
    /*<<<< transporters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    

    /*--------------------------------------------------------------------- */
    /* UPDATE AMPARs (set hits to 0) */
    /*>>>> AMPAR, set ampar_hits to 0 */
        memset(ampar_hits, 0, amparNo*sizeof(int));
        /* for (i_ampar = 0; i_ampar < amparNo; i_ampar++)
            *(ampar_hits + i_ampar) = 0; */
    /*<<<< AMPAR, set ampar_hits to 0 */
    
    /* UPDATE AMPARs (positions) */
    /* */
    /*>>>> AMPAR diffusion */
    if (D_ampar_PSD != 0  ||  D_ampar_outside != 0) {
        
        /* NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE
         * the following code moves an AMPAR through the PSD, cleft, ES,
         * but this is only valid as long as the mean square displacement dt_ampar
         * is small compared to the PSD radius (reason: faster performance). 
         * This should be true for moderately unphysiological values of D_ampar.
         */
        
        for (i_ampar = 0; i_ampar < amparNo; i_ampar++) {
            /* v = dt_ampar*randn(1, 3); */
            if (*(ampar_tag + i_ampar) == amparState_InCleft_OnPSD) {
                v[0] = dt_ampar_PSD     *gaussrand();
                v[1] = dt_ampar_PSD     *gaussrand();
            } else {
                v[0] = dt_ampar_outside *gaussrand();
                v[1] = dt_ampar_outside *gaussrand();
            }
            
            /* v(3) = 0; /* don't change 3d coordinate---diffusion is 2D in PSD */

            /* X = last valid position
               Y = suggested new position */
            /* X = ampar.traj(:, i_ampar); */
                X1 = *(ampar_traj + 0 + i_ampar*3);
                X2 = *(ampar_traj + 1 + i_ampar*3);
                X3 = *(ampar_traj + 2 + i_ampar*3);
                
            /* Y = X + v; */
                Y1 = X1 + v[0];
                Y2 = X2 + v[1];
                Y3 = X3;

            checking = 1;
            while (checking == 1) {
                checking = 0;

                /* where the transition into another area (PSD, outside-cleft, ES)
                 * is made, there the new (proposed) position vector (v resp Y)
                 * is determined.
                 */

                if (*(ampar_tag + i_ampar) == amparState_InCleft_OnPSD) {
                    
                    /**** transit to outside of PSD ****/
                    R = R_PSD; /* use here R_PSD / R_Cleft / R_Reservoir */
                    if (Y1*Y1 + Y2*Y2 > R*R) {
                        checking = 1;
                        
                        *(ampar_tag + i_ampar) = amparState_InCleft_OutsidePSD;
                        
                        theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                        + sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                        + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                        - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                        /* Z is the point where the segment X to Y intersects the */
                        /* PSD boundary */
                        Z1 = theta*Y1 + (1-theta)*X1;
                        Z2 = theta*Y2 + (1-theta)*X2;

                        /* now compare next point with last point inside, that's Z */
                        X1 = Z1;
                        X2 = Z2;
                        
                        /* reflect at PSD boundary with probability of P_Reflect_Inside,
                           otherwise leave to outside of PSD */
                        if (rand() <= P_Reflect_Inside) {
                            distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2;
                            Y1 = Y1 - 2*distTimesR*Z1/(R*R);
                            Y2 = Y2 - 2*distTimesR*Z2/(R*R);
                            
                            *(ampar_tag + i_ampar) = amparState_InCleft_OnPSD;
                        } /* reflect at PSD */
                        else { 
                            /* transit to PSD*/
                            v[0] = Y1 - X1;
                            v[1] = Y2 - X2;
                            /* rescale from InCleft_OnPSD diffusion to
                               InCleft_OutsidePSD diffusion */
                            Y1 = X1 + v[0]*dt_ampar_outside/dt_ampar_PSD;
                            Y2 = X2 + v[1]*dt_ampar_outside/dt_ampar_PSD;
                        }
                    }
                } /* AMPAR on PSD */
                
                if (*(ampar_tag + i_ampar) == amparState_InCleft_OutsidePSD) {
                    
                    /**** transit to ES reservoir ****/
                    R = R_Cleft; /* use here R_PSD / R_Cleft / R_Reservoir */
                    if (Y1*Y1 + Y2*Y2 > R*R) {
                        checking = 1;
                        
                        *(ampar_tag + i_ampar) = amparState_InES;
                        
                        theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                        + sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                        + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                        - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                        /* Z is the point where the segment X to Y intersects the */
                        /* PSD boundary */
                        Z1 = theta*Y1 + (1-theta)*X1;
                        Z2 = theta*Y2 + (1-theta)*X2;
                        
                        /* distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2; */
                        X1 = Z1;
                        X2 = Z2;
                    }
                    
                    /**** transit to PSD ****/
                    /* find where new (reflected) position goes back into PSD */
                    R = R_PSD; /* use here R_PSD / R_Cleft / R_Reservoir */
                    if (Y1*Y1 + Y2*Y2 < R*R) {
                        checking = 1;

                        *(ampar_tag + i_ampar) = amparState_InCleft_OnPSD;
                        
                        theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                        - sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                        + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                        - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                        
                        Z1 = theta*Y1 + (1-theta)*X1;
                        Z2 = theta*Y2 + (1-theta)*X2;
                        
                        X1 = Z1;
                        X2 = Z2;
                        
                        /* reflect at PSD boundary with probability of P_Reflect_Outside,
                           otherwise enter PSD */
                        if (rand() <= P_Reflect_Outside) {
                            distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2;
                            Y1 = Y1 - 2*distTimesR*Z1/(R*R);
                            Y2 = Y2 - 2*distTimesR*Z2/(R*R);
                            
                            *(ampar_tag + i_ampar) = amparState_InCleft_OutsidePSD;
                        } /* if reflect ... */
                        else { 
                            /* transit to PSD*/
                            v[0] = Y1 - X1;
                            v[1] = Y2 - X2;
                            /* rescale from InCleft_OutsidePSD diffusion to
                               InCleft_OnPSD diffusion */
                            Y1 = X1 + v[0]*dt_ampar_PSD/dt_ampar_outside;
                            Y2 = X2 + v[1]*dt_ampar_PSD/dt_ampar_outside;
                        }
                    }
                    
                } /* AMPAR outside PSD, inside cleft */
                
                if (*(ampar_tag + i_ampar) == amparState_InES) {
                    
                    /**** reflect at outer reservoir boundary ****/
                    R = R_Reservoir; /* use here R_PSD / R_Cleft / R_Reservoir */
                    if (Y1*Y1 + Y2*Y2 > R*R) {
                        checking = 1;
                        
                        theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                        + sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                        + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                        - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                        /* Z is the point where the segment X to Y intersects the */
                        /* PSD boundary */
                        Z1 = theta*Y1 + (1-theta)*X1;
                        Z2 = theta*Y2 + (1-theta)*X2;
                        
                        distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2;
                        /* reflect Y on the plane tangential in Z */
                        Y1 = Y1 - 2*distTimesR*Z1/(R*R);
                        Y2 = Y2 - 2*distTimesR*Z2/(R*R);
                        /* now compare next point with last point inside, that's Z */
                        X1 = Z1;
                        X2 = Z2;
                    }
                 
                    /**** transit to cleft ****/
                    /* find where new (reflected) position goes back into cleft */
                    R = R_Cleft; /* use here R_PSD / R_Cleft / R_Reservoir */
                    if (Y1*Y1 + Y2*Y2 < R*R) {
                        checking = 1;

                        *(ampar_tag + i_ampar) = amparState_InCleft_OutsidePSD;
                        
                        theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                        - sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                        + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                        - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                        
                        Z1 = theta*Y1 + (1-theta)*X1;
                        Z2 = theta*Y2 + (1-theta)*X2;
                        
                        /* distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2; */
                        X1 = Z1;
                        X2 = Z2;
                    }
                    
                } /* AMPAR in ES, not in cleft */
                
                
            } /* while checking */

            *(ampar_traj + 0 + i_ampar*3) = Y1;
            *(ampar_traj + 1 + i_ampar*3) = Y2;
            *(ampar_traj + 2 + i_ampar*3) = Y3;
        } /* for i_ampar */
    } /* if diffusion_AMPAR = ON ... */

    /*<<<< AMPAR diffusion */

    
    /*---------------------------------------------------------------------- */
    /* UPDATE glu's */
    /* */
    /*>>>> glu diffusion */

/*    if (*(gluStates + timestep + (gluState_FreeCleft-1)*timeSteps) > 0  ||
        *(gluStates + timestep + (gluState_UnbindingFromTransp-1)*timeSteps) > 0  ||
        *(gluStates + timestep + (gluState_UnbindingFromAMPAR-1)*timeSteps) > 0  ||
        *(gluStates + timestep + (gluState_FreeES-1)*timeSteps) > 0) {
 */
    if (1) {
        for (i_glu = 0; i_glu < gluNo; i_glu++) {
            int gluHasJustUnbound = 0;

            /* get new displacement vector v upon unbinding from a transporter */
            if (*(glu_state + i_glu) == gluState_UnbindingFromTransp) {
                double  R = *(insertDistrTable + (rand() & insertDistrTableLen));
                double  transpX = *(glu_traj + 0 + i_glu*3); /* this is the transporter X coord to which the glu is bound */
                double  transpY = *(glu_traj + 1 + i_glu*3); /* Y coord */
                double  l = R_ES; /* R_ES = sqrt(transpX*transpX + transpY*transpY); */
                v[0] = -R* transpX / l; /* release inwards */
                v[1] = -R* transpY / l;
                v[2] = dt_glu*gaussrand();

                *(glu_state + i_glu) = gluState_FreeES;

                gluHasJustUnbound = 1;
            }

            /* get new displacement vector v upon unbinding from an AMPAR */
            if (*(glu_state + i_glu) == gluState_UnbindingFromAMPAR) {
                double  R = *(insertDistrTable + (rand() & insertDistrTableLen));
                double  transpX = *(glu_traj + 0 + i_glu*3); /* this is the transporter X coord to which the glu is bound */
                double  transpY = *(glu_traj + 1 + i_glu*3); /* Y coord */
                
                *(glu_state + i_glu) = gluState_FreeCleft;

                v[0] = dt_glu*gaussrand();
                v[1] = dt_glu*gaussrand();
                v[2] = R;
                
                *(glu_state + i_glu) = gluState_FreeCleft;

                gluHasJustUnbound = 1;
            }
                
            if (*(glu_state + i_glu) == gluState_FreeCleft  ||  *(glu_state + i_glu) == gluState_FreeES) {
                if (gluHasJustUnbound == 0) { /* get new displacement vector v*/
                    v[0] = dt_glu*gaussrand();
                    v[1] = dt_glu*gaussrand();
                    v[2] = dt_glu*gaussrand();
                }
                /* X = glu.traj(:, i_glu); */
                X1 = *(glu_traj + 0 + i_glu*3);
                X2 = *(glu_traj + 1 + i_glu*3);
                X3 = *(glu_traj + 2 + i_glu*3);

                /* Y = X + v; */
                Y1 = X1 + v[0];
                Y2 = X2 + v[1];
                Y3 = X3 + v[2];
                
                /* in the following: 
                 *   X = last valid position
                 *   Y = suggested new position
                 */

                /* The variable `checking ' controls whether the particle position
                 * is updated and a collision check is performed */
                checking = 1; /* update and collision-detect */

                /* position update and collision detection */
                while (checking == 1) {

                    /*>>>> glu in cleft, not in extrasynaptic space */
                    if (*(glu_state + i_glu) == gluState_FreeCleft) {

                        checking = 0;

                        /*>>>> absorb on leaving the cleft
                          if (Y1*Y1 + Y2*Y2 > R_Cleft*R_Cleft) {
                              *(glu_state + i_glu) = gluState_Absorbed;
                              *(glu_exitTime + i_glu) = timestep*timeStepSize;
                              checking = 0;
                              break;
                          }
                         */
                        /*<<<< absorb on leaving the cleft*/
                        
                        /* XOR */
                        
                        /*>>>> transit to extrasynaptic space on hitting
                         * the lateral cleft boundary
                         * If Y will leave the cleft, we will need further checking
                         * in the extrasynaptic space */
                        if (Y1*Y1 + Y2*Y2 > R_Cleft*R_Cleft) {
                            *(glu_state + i_glu) = gluState_FreeES; /* update position flag already here */
                            checking = 1;
                        }
                        /*<<<< transit to extrasynaptic space */


                        /* reflect at cylinder top and bottom */
                        /* X is last valid point in the space */
                        /* this will be updated, upon reflection, to the point of reflection */
                        Z1 = X1;
                        Z2 = X2;
                        Z3 = X3;
                        while (Y3 > upperLid  ||  Y3 < lowerLid) {
                            if (Y3 > upperLid) {
                                /* Z = point of hitting the PLANE through the upper lid*/
                                Z3 = upperLid;
                                Z1 = X1 + (Z3-X3)*(Y1-X1)/(Y3-X3);
                                Z2 = X2 + (Z3-X3)*(Y2-X2)/(Y3-X3);
                                if (Z1*Z1 + Z2*Z2 < R_Cleft*R_Cleft) {
                                    Y3 = 2*upperLid-Y3; /* reflect */
                                    X1 = Z1;
                                    X2 = Z2;
                                    X3 = Z3;
                                } else { /* jump into ES */
                                    break; /* new position flag has been set already */
                                }
                            }
                            if (Y3 < lowerLid) {
                                /* Z = point of hitting the PLANE through the lower lid*/
                                Z3 = lowerLid;
                                Z1 = X1 + (Z3-X3)*(Y1-X1)/(Y3-X3);
                                Z2 = X2 + (Z3-X3)*(Y2-X2)/(Y3-X3);
                                
                                if (Z1*Z1 + Z2*Z2 < R_Cleft*R_Cleft) {
                                    Y3 = 2*lowerLid-Y3; /* reflect */
                                    X1 = Z1;
                                    X2 = Z2;
                                    X3 = Z3;
                                    
                                /*>>>> AMPAR state transition on the event of hitting ++++++++++++++++++++++*/
                                    for (i_ampar = 0; i_ampar < amparNo; i_ampar++) {
                                        double  P1 = *(ampar_traj + 0 + i_ampar*3); /* AMPAR position */
                                        double  P2 = *(ampar_traj + 1 + i_ampar*3);
                                        /* take only AMPARs inside the clect */
                                        if (P1*P1 + P2*P2 < R_Cleft*R_Cleft) {
                                            double  A1 = Z1 - P1;
                                            double  A2 = Z2 - P2;
                                            
                                            if (A1*A1 + A2*A2 < R_BindingToAMPAR*R_BindingToAMPAR) {
                                                (*(ampar_hits + i_ampar)) ++;
                                                
                                            /* we DO NOT BIND ligands in this version !!! */
                                            } /* if glu in R_BindingToAMPAR */
                                        }
                                
                                    } /* for i_ampar */
                                /*<<<< AMPAR state transition on the event of hitting ++++++++++++++++++++++*/
                                    
                                } else { /* jump into ES */
                                    break; /* new position flag has been set already */
                                }
                            }
                        } /* while reflect at top/bottom plane */

#ifdef ABSORB_AT_CLEFT_BOUNDARY
                        if (Y1*Y1 + Y2*Y2 > R_Cleft*R_Cleft) {
                            checking = 0;
                            checkingES = 0;
                            *(glu_state + i_glu) = gluState_Absorbed;
                        }
#endif
                        
                    } /* if glu.state == gluState_FreeCleft */
                    /*<<<< glu in cleft, not in extrasynaptic space */


                    /*>>>> glu in extrasynaptic space, not in cleft */
                    if (*(glu_state + i_glu) == gluState_FreeES) {
                        checking = 0;
                        checkingES = 1;
                        while (checkingES == 1) {
                            checkingES = 0;

                            /*>>>> absorb on hitting the ES lateral boundary
                             * if (Y(1)^2 + Y(2)^2 > R_ES^2) {
                             *     glu.state(i_glu) = gluState_Absorbed;
                             *     glu.exitTime(i_glu) = timestep*timeStepSize;
                             *     checking = 0;
                             * }
                             *<<<< absorb on hitting the ES lateral boundary
                             * OR
                             *>>>> reflection on hitting outer ES lateral boundary */
                            R = R_ES; /* dont remove, R is used in formulas below */
                            if (Y1*Y1 + Y2*Y2 > R*R) {
                                
#ifdef ABSORB_AT_GLIAL_SHEATH
                                checking = 0;
                                checkingES = 0;
                                *(glu_state + i_glu) = gluState_Absorbed;
                                break;
#endif
                                
                                checkingES = 1;

                                theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                                        + sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                                             + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                                             - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                                /* Z is the point where the segment X to Y intersects the */
                                /* ES lateral boundary */
                                Z1 = theta*Y1 + (1-theta)*X1;
                                Z2 = theta*Y2 + (1-theta)*X2;
                                /* division by zero cannot occur, as Y outside and X inside ES */
                                Z3 = X3 + sqrt( ((Z1-X1)*(Z1-X1)+(Z2-X2)*(Z2-X2))/((Y1-X1)*(Y1-X1)+(Y2-X2)*(Y2-X2)) ) *(Y3-X3);

                                distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2;
                                /* reflect Y on the plane tangential in Z */
                                Y1 = Y1 - 2*distTimesR*Z1/(R*R);
                                Y2 = Y2 - 2*distTimesR*Z2/(R*R);
                                /* set `last valid point', ie X, to point of hitting, ie Z */
                                X1 = Z1;
                                X2 = Z2;
                                X3 = Z3;
                                
                                /*>>>> binding to transporter ++++++++++++++++++++++++++++++++++++++++++++++++*/
                                /* Test hitting on transporter */
                                if (transpNo > 0) {
                                    /* find the right transporter ring(s) */
                                    int i_ver_U = (int)floor( (Z3+R_BindingToTransp - lowerLidES)/transp2transpDist );
                                    int i_ver_D = (int)floor( (Z3-R_BindingToTransp - lowerLidES)/transp2transpDist );
                                    
                                    for (i_ver = i_ver_D; i_ver < i_ver_U; i_ver++) {
                                        if (i_ver >= 0  &&  i_ver < transpVer) {
                                            for (i_hor = 0; i_hor < transpHor; i_hor++) {
                                                
                                                
                                                if (rand() <= P_bindingToTransp) { /* test whether hit counts */                                                    
                                                    double P1, P2, P3;
                                                    
                                                    if (*(transp_state + i_hor + i_ver*transpHor) != transpState_unbound) {
                                                        continue;
                                                    }
                                                    
                                                /* P = [Z1, Z2, Z3] - transp.pos(:, i_hor, i_ver); */
                                                    P1 = Z1 - *(transp_pos + 0 + (i_hor +  i_ver*transpHor)*3);
                                                    P2 = Z2 - *(transp_pos + 1 + (i_hor +  i_ver*transpHor)*3);
                                                    P3 = Z3 - *(transp_pos + 2 + (i_hor +  i_ver*transpHor)*3);
                                                    
                                                    if (P1*P1 + P2*P2 + P3*P3 < R_BindingToTransp*R_BindingToTransp) {
                                                        *(transp_ligand + i_hor + i_ver*transpHor) = i_glu;

                                                        /* find transition to the next state */
                                                        
                                                        if ( rand() < RAND_MAX*transpRate_pump/(transpRate_pump + transpRate_dissoc) ) {
                                                            /* (stay associated and) pump in */
                                                            *(transp_state  + i_hor + i_ver*transpHor) = transpState_bound_nextState_pumpedIn;
                                                                /* deterministic transporters */
                                                                /* *(transp_T      + i_hor + i_ver*transpHor) = transp_BindingTime; */
                                                            /* stochastic transporters */
                                                            *(transp_T      + i_hor + i_ver*transpHor) = -log((double)rand()/RAND_MAX) / transpRate_pump;
                                                        } else {
                                                            /* dissociate */
                                                            *(transp_state  + i_hor + i_ver*transpHor) = transpState_bound_nextState_unbound;
                                                                /* deterministic transporters */
                                                                /* *(transp_T      + i_hor + i_ver*transpHor) = transp_UnbindingTime; */
                                                            /* stochastic transporters */
                                                            *(transp_T      + i_hor + i_ver*transpHor) = -log((double)rand()/RAND_MAX) / transpRate_dissoc;
                                                        }
                                                            
                                                        /* binding: update glu */
                                                        *(glu_state + i_glu) = gluState_BoundToTransp;
                                                        /*glu.traj(:, i_glu) = transp.pos(:, i_hor, i_ver); */
                                                        Y1 = *(transp_pos + 0 + (i_hor + i_ver*transpHor)*3);
                                                        Y2 = *(transp_pos + 1 + (i_hor + i_ver*transpHor)*3);
                                                        Y3 = *(transp_pos + 2 + (i_hor + i_ver*transpHor)*3);
/*
breakit = 1;
 */
                                                        checking = 0; /* don't update glu[i_glu] any further */
                                                        checkingES = 0; /* don't update glu[i_glu] any further */
                                                        break; /* don't continue testing for binding */
                                                    }
                                                } /* if hit counts, ie, rand() < P_binding... */
                                            } /* for all transporters in the ring */
                                            
                                                /* if glu was bound to transporter, break `for transporter rings' loop */
                                            if (*(glu_state + i_glu) == gluState_BoundToTransp)
                                                break;
                                            
                                        } /* if glu in right transporter ring */
                                    } /* for some transporter rings... */
                                    
                                } /* if transpNo > 0 */
                                /*<<<< binding to transporter ++++++++++++++++++++++++++++++++++++++++++++++++*/
                                
                            } /* if Y outside ES lateral boundary */
                            /*<<<< reflection on hitting outer ES lateral boundary */

                            /*>>>> reflection on hitting inner ES lateral boundary */
                            R = R_Cleft; /* dont remove, R is used in formulas below */
                            if (Y1*Y1 + Y2*Y2 < R*R) {
                                checkingES = 1;

                                theta = (-X1*Y1-X2*Y2+X1*X1+X2*X2
                                        - sqrt(-2*X1*Y1*R*R-2*X2*Y2*R*R-X1*X1*Y2*Y2
                                             + 2*X1*Y1*X2*Y2+X1*X1*R*R+Y1*Y1*R*R
                                             - X2*X2*Y1*Y1+X2*X2*R*R+Y2*Y2*R*R))
                                        /(X1*X1-2*X1*Y1+Y1*Y1+X2*X2-2*X2*Y2+Y2*Y2);
                                
                                Z1 = theta*Y1 + (1-theta)*X1;
                                Z2 = theta*Y2 + (1-theta)*X2;
                                /* division by zero cannot occur, as Y outside ES and X inside ES */
                                Z3 = X3 + sqrt( ((Z1-X1)*(Z1-X1)+(Z2-X2)*(Z2-X2))/((Y1-X1)*(Y1-X1)+(Y2-X2)*(Y2-X2)) ) *(Y3-X3);
                                distTimesR = (Y1-Z1)*Z1 + (Y2-Z2)*Z2;
                                X1 = Z1;
                                X2 = Z2;
                                X3 = Z3;

                                /* Test if point Z of reflection lies in entry to cleft:
                                 * If transition to cleft, continue checking in cleft */
                                if (Z3 <= upperLid  &&  Z3 >= lowerLid) {
                                    *(glu_state + i_glu) = gluState_FreeCleft; /* update position flag */
                                    checking = 1;
                                    break; /* leave `while checkingES...', but continue `while checking...' loop */
                                }

                                /* reflect at inner ES boundary */
                                Y1 = Y1 - 2*distTimesR*Z1/(R*R);
                                Y2 = Y2 - 2*distTimesR*Z2/(R*R);
                            } /* if */
                            /*<<<< reflection on hitting inner ES lateral boundary */

                        } /* while checkingES */

                        /*>>>> absorb at ES cylinder top and bottom */
                        if (Y3 >= upperLidES  ||  Y3 <= lowerLidES) {
                            *(glu_state + i_glu) = gluState_Absorbed;
                            *(glu_exitTime + i_glu) = timestep*timeStepSize;
                            checking = 0;
                            checkingES = 0;
                            break; /* leave `while checking...' loop */
                        }
                        /*<<<< */

                    } /* if glu.state == gluState_FreeES */
                    /*<<<< glu in extrasynaptic space, not in cleft */

                } /* while checking */

                /* report glu in gluDistrES */

                if (*(glu_state + i_glu) == gluState_FreeES) {                    
                    int esBin = (int)floor( (Y3 - lowerLidES)/esBinSize );
                    if (esBin >= 0  &&  esBin < esBinNo) {
                        (*(gluDistrES + timestep + esBin*timeSteps)) ++;
                    }
                }

                
                *(glu_traj + 0 + i_glu*3) = Y1; /* update position */
                *(glu_traj + 1 + i_glu*3) = Y2; /* update position */
                *(glu_traj + 2 + i_glu*3) = Y3; /* update position */
            } /* if glu.state = Free and alive ... */

        } /* for i_glu */
    } /* if no glus in cleft and ES */
    /*<<<< glu diffusion */

    
    /*---------------------------------------------------------------------- */
    /* UPDATE AMPAR states */
    /* */
    /*>>>> AMPAR state transition >>>>>>>>>>>>>>>>>>>> */
    {
        double  P[MAX_AMPARSTATES_NO], Po[MAX_AMPARSTATES_NO], K[MAX_AMPARSTATES_NO], S[MAX_AMPARSTATES_NO];
        double  ref, t, h, H, to, Ho, errorSq;
        int     i, j, k, l;
        double  concentration;
        double  X1, X2, X3, r;
                                    
        double  s[MAX_AMPARSTATES_NO][6];
        
        int     kineticScheme;

        
    /* update all AMPARs */
    for (i_ampar = 0; i_ampar < amparNo; i_ampar++) {
        
        kineticScheme = *(ampar_kinetics + i_ampar) - 1; /* matlab indexing */
        

        /* r = distance from release (at 0)*/
                X1 = *(ampar_traj + 0 + i_ampar*3);
                X2 = *(ampar_traj + 1 + i_ampar*3);
                X3 = *(ampar_traj + 2 + i_ampar*3);
        r = sqrt(X1*X1 + X2*X2);

        memcpy(P, (ampar_stateDistr + i_ampar *MAX_AMPARSTATES_NO), MAX_AMPARSTATES_NO*sizeof(double));        

        t = timeStepSize*timestep; /* actual time */
                
        if (*(ampar_hits + i_ampar) == 0) {
    
            /* matrix multiplication */
            for (i_states = 0; i_states < MAX_AMPARSTATES_NO; i_states++) {
    
                P[i_states] = 0;
                for (i = 0; i < MAX_AMPARSTATES_NO; i++)
                    P[i_states] += P_transition[kineticScheme][i][i_states] * (*(ampar_stateDistr + i + i_ampar *MAX_AMPARSTATES_NO));
            }
                                    
        } else { /* use RK solver */

            concentration = (double)*(ampar_hits + i_ampar) * hitsToConcentration;

            /*>>>> advance time of timeStepSize */
            /*Runge-Kutta solver */
            ref = 1; /* no refinement */
            h = timeStepSize*ref;
            H = 0; /* portion of time step already calculated */
            
            while (H < 1) {
                memcpy(Po, P, MAX_AMPARSTATES_NO*sizeof(double)); /* Po = P; */
                to = t;
                Ho = H;
                
                /* get new intermediate values s[] */
                for (i = 0; i < 6; i++) {
                    memcpy(S, P, MAX_AMPARSTATES_NO*sizeof(double)); /* S = P; */
                    for (j = 0; j < i-1; j++)
                        for (k = 0; k < MAX_AMPARSTATES_NO; k++)
                            S[k] += RK_B[i][j]*s[k][j];
                    
                 /* get time-dependent, local concentration here */
                 /* current time at this point: T = t + RK_A[i]*h; */

/*
                * update concentration-sensitive transitions *
                for (i_1, j_1 = 0; i < MAX_AMPARSTATES_NO; i++) {
                    Q_transition[i_1][j_1] = <<RATE>> *concentration;
                }
                * update diagonal *
                for (i = 0; i < MAX_AMPARSTATES_NO; i++) {
                    Q_transition[i][i] = <<old_diagonal>> + <<updated_transitions>> ;
                }
*/                    
                    
                    
                if (kineticScheme == kinetics_a1 - 1) {                    
                    /* update concentration-sensitive transitions */
                    Q_transition[kineticScheme][amparState_C0_a1-1][amparState_C1_a1-1] = k1_a1 *concentration;
                    Q_transition[kineticScheme][amparState_C1_a1-1][amparState_C2_a1-1] = k2_a1 *concentration;
                    Q_transition[kineticScheme][amparState_D3_a1-1][amparState_D4_a1-1] = k5_a1 *concentration;
                    
                    /* update diagonal */
                    Q_transition[kineticScheme][amparState_C0_a1-1][amparState_C0_a1-1] =            - Q_transition[kineticScheme][amparState_C0_a1-1][amparState_C1_a1-1];
                    Q_transition[kineticScheme][amparState_C1_a1-1][amparState_C1_a1-1] = Q_C1_C1_a1 - Q_transition[kineticScheme][amparState_C1_a1-1][amparState_C2_a1-1];
                    Q_transition[kineticScheme][amparState_D3_a1-1][amparState_D3_a1-1] = Q_D3_D3_a1 - Q_transition[kineticScheme][amparState_D3_a1-1][amparState_D4_a1-1];
                } else {
                    /* update concentration-sensitive transitions */
                    Q_transition[kineticScheme][amparState_C0_a4-1][amparState_C1_a4-1] = k1_a4 *concentration;
                    Q_transition[kineticScheme][amparState_C1_a4-1][amparState_C2_a4-1] = k2_a4 *concentration;
                    Q_transition[kineticScheme][amparState_D3_a4-1][amparState_D4_a4-1] = k5_a4 *concentration;
                    
                    /* update diagonal */
                    Q_transition[kineticScheme][amparState_C0_a4-1][amparState_C0_a4-1] =            - Q_transition[kineticScheme][amparState_C0_a4-1][amparState_C1_a4-1];
                    Q_transition[kineticScheme][amparState_C1_a4-1][amparState_C1_a4-1] = Q_C1_C1_a4 - Q_transition[kineticScheme][amparState_C1_a4-1][amparState_C2_a4-1];
                    Q_transition[kineticScheme][amparState_D3_a4-1][amparState_D3_a4-1] = Q_D3_D3_a4 - Q_transition[kineticScheme][amparState_D3_a4-1][amparState_D4_a4-1];
                }
                    
                    /* S = h*S */
                    for (k = 0; k < MAX_AMPARSTATES_NO; k++)
                        S[k] *= h;
                    
                    /* s(:, i) = h * transpose(Q) * S; */
                    for (k = 0; k < MAX_AMPARSTATES_NO; k++) {
                        s[k][i] = 0;
                        for (l = 0; l < MAX_AMPARSTATES_NO; l++)
                            s[k][i] += Q_transition[kineticScheme][l][k] * S[l];
                    } /* for k */
                    
                } /* for i */
                
                /* get corrector K and predictor P */
                memcpy(K, P, MAX_AMPARSTATES_NO*sizeof(double)); /* K = P */
                for (i = 0; i < 6; i++)
                    for (k = 0; k < MAX_AMPARSTATES_NO; k++) {
                        K[k] += RK_C[0][i] * s[k][i];
                        P[k] += RK_C[1][i] * s[k][i];
                    }
                
                /* get norm of error */
                errorSq = 0;
                for (k = 0; k < MAX_AMPARSTATES_NO; k++)
                    errorSq += (P[k] - K[k])*(P[k] - K[k]); /* (P - K)^2 */
                
                /* adapt time step according to error */
                if (errorSq > RK_ErrorTol*RK_ErrorTol) {
                    /* P = Po; discard newly calculated data point */
                    memcpy(P, Po, MAX_AMPARSTATES_NO*sizeof(double));
                    
                    if (ref <= RK_MaxRef) {
                        t = t + h;
                        H = H + ref;
                    } else {
                        ref /= 2;
                        h = timeStepSize*ref;
                    }
                } else {
                    t = t + h;
                    H = H + ref;
                    if (ref < 1) {
                        ref *= 2;
                        h = timeStepSize*ref;
                    }
                }
                
            } /* while (H < 1) */
            
            /* interpolate endpoint
               `to' is last point inside timeStepSize interval <
               T := timeStepSize*(timestep+1) is endpoint of timeStepSize interval <
               `ta' is farther time point the newly calculated `Pa' corresponds to */
            if (H > 1) {
                double lambda = ( (timeStepSize*(timestep+1)) -to)/(t-to);
                for (k = 0; k < MAX_AMPARSTATES_NO; k++)
                    P[k] = Po[k] + (P[k]-Po[k])*lambda;
            }
            /*<<<< advance time of timeStepSize */
        } /* use RK solver */

        /* update state distribution of AMPAR number i_ampar */
        memcpy((ampar_stateDistr + i_ampar *MAX_AMPARSTATES_NO), P, MAX_AMPARSTATES_NO*sizeof(double));        
        
        
    } /* for i_ampar */
                    

    }
    /*<<<< AMPAR state transition <<<<<<<<<<<<<<<<<<<< */

    
    
    
    /* UPDATE timestep */
    /* */
    /* increase timestep at the end of loop */
    timestep ++;
    
    
    /*#####################################################################*/
    {
        /* double firstDec; */
        /* percent calculated */
        /*if (mod(timestep, round(timeSteps/100*1)) == 0) */
        /*firstDec = floor(((timestep/timeSteps)*100-floor(timestep/timeSteps*100))*10);*/
        
        mexCallMATLAB(0, NULL, 0, NULL, "printLoopStateInfo");
    /*        mexPrintf("\b\b\b\b\b%3.0f.%1.0f", floor((timestep/timeSteps)*100), firstDec * (firstDec < 10)); */
        /*end */
    }

    
} /* while timestep ... */

/*  matrix access:
    a(i1, i2, i3) = *(a_ptr + i1 + (i2 + (i3)*i2_size)*i1_size)
 */

/*  {
    double  double_dummy;
    var_ptr = getMATLABvar("timestep_internal");
    double_dummy = timestep + 1;
    mxSetPr(var_ptr, &double_dummy);
    setMATLABvar_ptr("timestep_internal", var_ptr);
    }
*/
    /* update book-keeping variables */
    setMATLABvar_ptr("amparTags",    amparTags_ptr);
    setMATLABvar_ptr("amparTagsVar", amparTagsVar_ptr);
    setMATLABvar_ptr("amparStates",  amparStates_ptr);
    setMATLABvar_ptr("amparMean",    amparMean_ptr);
    setMATLABvar_ptr("amparVarS",    amparVarS_ptr);
    setMATLABvar_ptr("gluStates",    gluStates_ptr);
    setMATLABvar_ptr("gluDistrES",   gluDistrES_ptr);

/*
    mxDestroyArray(amparTags_ptr);
    mxDestroyArray(amparTagsVar_ptr);
    mxDestroyArray(amparStates_ptr);
    mxDestroyArray(amparMean_ptr);
    mxDestroyArray(amparVarS_ptr);
    mxDestroyArray(gluStates_ptr);
    mxDestroyArray(gluDistrES_ptr);
*/    
} /* mexFunction */
