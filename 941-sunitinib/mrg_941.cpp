$PROB Sunitinib, Yu et al, Br J Clin Pharmacol, 2014

$PARAM @annotated
K12   : 0.34 : Absorption rate constant 12 (h-1)
TVCLS : 35.7 : Typical Value Clearance Sunitinib (L/h)
TVVCS : 1360 : Typical Value Central Volume Sunitinib (L)
TVCLM : 17.1 : Typical Value Clearance Metabolite (L/h)
TVVCM : 635  : Typical Value Central Volume Metabolite (L)
TVQ   : 20.1 : Typical Value Intercompartimental Clearance Metabolite (L/h)
TVVPM : 388  : Typical Value Peripheral volume Metabolite (L)

ETA1 : 0 : Central Volume Sunitinib(L)
ETA2 : 0 : Central Volume Metabolite (L)
ETA3 : 0 : Clearance metabolite (L/h)
ETA4 : 0 : Clearance sunitinib (L/h)

$PARAM @annotated @covariates
WT : 70 : Body weight (kg)

$CMT @annotated
DEPOT : Depot [ADM]
CENTRALS : Central sunitinib [OBS]
CENTRALM : Central metabolite [OBS]
PERIPHM : Peripheral metabolite

$OMEGA @block 
0.104976 // VC suni
0.016892 0.335241 // VC metab
0.000000 0.026738 0.177241 // CL metab
0.000000 0.000000 0.010795 0.114921 // CL suni

$SIGMA
0.043 0 // prop suni
0.023 0 // prop metab

$MAIN
double ASCL = pow((WT/70), 0.75) ;
double ASV  = WT/70;
double FM   = 0.21 ;
double QH   = 80*ASCL ;

double CLS = TVCLS * ASCL * exp(ETA(4) + ETA4) ; 
double VCS = TVVCS * ASV  * exp(ETA(1) + ETA1) ;

double CLM = TVCLM * ASCL * exp(ETA(3) + ETA3) ;
double VCM = TVVCM * ASV  * exp(ETA(2) + ETA2) ; 

double VPM = TVVPM * ASV ;
double Q = TVQ * ASCL ;
double K34 = Q / VCM  ;
double K43 = Q / VPM  ;

$ODE
double CLIV = ((K12 * DEPOT) + (QH * CENTRALS / VCS)) / (QH + CLS) ;

dxdt_DEPOT    = - K12 * DEPOT ;
dxdt_CENTRALS = QH * CLIV  - QH * CENTRALS / VCS ;
dxdt_CENTRALM = FM * CLS * CLIV - CLM * CENTRALM / VCM - K34 * CENTRALM + K43 * PERIPHM ;
dxdt_PERIPHM  = K34 * CENTRALM - K43 * PERIPHM ;

$ERROR
double PAR = 1000 * (CENTRALS / VCS) * (1 + EPS(1)) + EPS(2) ; 
double MET = 1000 * (CENTRALM / VCM) * (1 + EPS(3)) + EPS(4) ;
double DV = PAR ;
if(self.cmt == 3) DV = MET ;

$CAPTURE DV PAR MET