$PROB Methotrexate, Faltaos et al, CCP, 2006

$PARAM @annotated
TVCL : 7.1 : Clearance (L/h)
TVVC : 25.1 : Central volume of distribution (L)
TVQ  : 0.15 : Intercompartmental clearance (L/h)
TVVP : 2.7 : Peripheral volume (L)

AGE_CL : -0.22 : Age on CL ()
SCR_CL : -0.43 : Serum creat on CL ()

ETA1 : 0 : CL
ETA2 : 0 : VC
ETA3 : 0 : Q
ETA4 : 0 : VP
ETA5 : 0 : IOVCL1
ETA6 : 0 : IOVCL2

$PARAM @annotated @covariates
AGE : 62 : Age (years)
SCR : 67 : Creatinine clearance (mL/min/1.73m2)
CYCLE : 1 : Cycle

$OMEGA @block
0.266256 // CL
0.002401 0.233289 // VC
$OMEGA @block
0.430336 // Q
0.084164 0.4096  // VP
$OMEGA
0.027225 // IOV CL 1
0.027225 // IOV CL 2

$SIGMA 
0.2116 // err prop
0.000225 //  err additive .15^2

$CMT @annotated
CENTRAL : Central () [ADM, OBS]
PERIPH : Periph ()

$TABLE
double DV = 2.2 * (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ; //amt mg, dv microM

$MAIN
double IOVCL = ETA(5) + ETA5 ;
if(CYCLE != 1) IOVCL = ETA(6) + ETA6 ;

double CL = TVCL * exp(ETA(1) + ETA1 + IOVCL) * pow((SCR / 67), SCR_CL) * pow((AGE / 62), AGE_CL) ; 
double VC = TVVC * exp(ETA(2) + ETA2 ) ;
double Q  = TVQ  * exp(ETA(3) + ETA3 ) ;
double VP = TVVP * exp(ETA(4) + ETA4 ) ;

double K20 = CL / VC ;
double K23 = Q / VC ;
double K32 = Q / VP ;

$ODE
dxdt_CENTRAL = - (K20 + K23) * CENTRAL + K32 * PERIPH ;
dxdt_PERIPH  = - K32 * PERIPH + K23 * CENTRAL ;  

$CAPTURE DV