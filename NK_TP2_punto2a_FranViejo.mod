%% TNT SGU - Rotemeberg prices

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

/// Endogenous variables-27total(8nuevas)
var c laC laW h w Rst ppiS ppi R cN pN p cT ppiT rer ppiN mcN yN m cW cC hW  mW tao ;
var gdp tb dst dstC d dC;

// Variables that hav associated an exogenous process
var z ppist em Rw yT ppiTar g;

//37var en total - 30 exogenas

// Disturbances variables (7)
varexo u_z u_ppist u_em u_Rw u_yT u_ppiTar u_g;

// Parameters 
parameters SIG VPHI CHI BET GAM VRHO NU XI;
parameters ALP EPS THET PHIN;
parameters ALP_PPI ALP_GDP BAR_PPI PHI_B BARD;
parameters RHO_z RHO_ppist RHO_em RHO_Rw RHO_yT RHO_ppiTar RHO_g;
parameters SIGMA_z SIGMA_ppist SIGMA_em SIGMA_Rw SIGMA_yT SIGMA_ppiTar SIGMA_g;
parameters DELTA PHID PHIG DEUDA;            //NUEVOS

parameters R_ss gdp_ss ppi_ss z_ss ppist_ss em_ss Rw_ss yT_ss ngdp_ss g_ss tao_ss; // in eqns and exon
parameters h_ss stb_ss rer_ss sm_ss syT_ss sg_ss d_ss ;                            // given in ss

%----------------------------------------------------------------
% 2. Calibration                
%----------------------------------------------------------------
// Calibarted parameteres and steady state values
% CALCULATED ENDOGENOUSLY
% ppist_ss, CHI, BARD, delF_ss, yT_ss, g_ss

BET          = 1/1.03;   // SGU
SIG          = 2;             // SGU
VPHI         = 0.5;           // SGU
VRHO         = 0.5;           // SGU
GAM          = 1-0.26;        // SGU

ALP        = 1;          // SGU     
EPS        = 6;             // SGU
THET       = 0.7;           // SGU

//Nuevos
DELTA=0.8;
PHID=0.3;
PHIG=0.1;


PHIN= (EPS-1)/((1-THET)*(1-THET*BET)/THET);

ALP_PPI      = 1.5;           // Taylor
ALP_GDP      = 0.5/4;         // Taylor
BAR_PPI      = 1;

RHO_z        = 0.9;
RHO_ppist    = 0.9;
RHO_em       = 0;
RHO_Rw       = 0.9;
RHO_yT       = 0.9;
RHO_ppiTar   = 0.9;
RHO_g       = 0.5;      //nuevo

SIGMA_z      = 0.01;
SIGMA_ppist  = 0.01;
SIGMA_em     = 0.01;
SIGMA_Rw     = 0.01;
SIGMA_yT     = 0.01;
SIGMA_ppiTar = 0.01;
SIGMA_g     = 0.01;

% NORMALIZATION EXOGENOUS VAR
em_ss    = 1;

z_ss =1;

Rw_ss    = 1.03;          // SGU

% EXON GIVEN
h_ss   = 1;

PHI_B   = 0.0000335;     // SGU 

BARD = 1;

yT_ss=0.4;

sm_ss=0.2;

XI=5;

sg_ss=0.1;

//NUevas variables de ss
DEUDA = 0;

%----------------------------------------------------------------
% 3. Model 
%----------------------------------------------------------------
model; 

// CPOs workers
exp(laW) = exp(cW)^(-SIG);                                             //E1
CHI*exp(hW)^VPHI = exp(laW)*exp(w);                                    //E2
exp(laW) = BET*exp(laW(+1))/exp(ppi(+1)) + NU*exp(mW)^(-XI);           //E3
exp(cW)+exp(mW)+exp(tao)=exp(w)*exp(hW)+exp(mW(-1))/exp(ppi);                   //E4

//CPOs capitalists
exp(laC) = exp(cC)^(-SIG);                                              //E5
exp(laC) = BET*exp(R)*exp(laC(+1))/exp(ppi(+1));                        //E6
exp(laC) = BET*exp(Rst)*exp(laC(+1))*exp(ppiS(+1))/exp(ppi(+1));        //E7

// Ecuaciones de agregación
exp(c)=DELTA*exp(cW)+(1-DELTA)*exp(cC);                                 //E8
exp(h)=DELTA*exp(hW);                                                   //E9
exp(m)=DELTA*exp(mW);                                                   //E10
dst=(1-DELTA)*dstC;                                                     //E11
d=(1-DELTA)*dC;                                                  //E12

//Ecuaciones de producción total
exp(cN) = GAM*(exp(pN)/exp(p))^(-VRHO)*exp(c);                          //E13
exp(cT) = (1-GAM)*(1/exp(p))^(-VRHO)*exp(c);                            //E14
exp(c) = (GAM^(1/VRHO)*exp(cN)^(1-1/VRHO) + (1-GAM)^(1/VRHO)*exp(cT)^(1-1/VRHO))^(VRHO/(VRHO-1));   //E15

//Monopolist
exp(yN)=exp(z)*exp(h)^(ALP);                                            //E16
exp(pN)*exp(rer)*exp(mcN)*exp(z)*ALP*exp(h)^(ALP-1)=exp(w);             //E17
EPS-1=exp(mcN)*(EPS)-PHIN*(exp(ppiN)-1)*exp(ppiN)+PHIN*BET*exp(laC(+1))/(exp(ppi(+1))*exp(laC))*exp(yN(+1))/exp(yN)*(exp(ppiN(+1))-1)*(exp(ppiN(+1)))^2;
//E18

//Eq tasas con premium
exp(Rst) = exp(Rw)*exp(PHI_B*(dst-BARD));                               //E19

//RER definition
exp(rer)=1/exp(p);                                                      //E20

//Mx clearing
exp(cN) = exp(yN)*(1-PHIN/2*(exp(ppiN)-1)^2)-exp(g);                    //E21

//Relative prices
exp(ppiT) = exp(rer)*exp(ppi)/exp(rer(-1));                             //E22
exp(ppiN) = exp(pN)*exp(ppiT)/exp(pN(-1));                              //E23
exp(rer)/exp(rer(-1))= exp(ppiS)*exp(ppist)/exp(ppi);                   //E24

// Mx completion international and aggregation
exp(Rst(-1))/exp(ppist)*dst(-1) = tb + dst;                             //E25
exp(gdp) = exp(c) + tb + exp(g);                                                 //E26
tb = exp(yT)-exp(cT);                                                   //E27

//Monetary Policy
// Taylor Rule
exp(R)= R_ss*exp(ppiTar)*(exp(ppi)/exp(ppiTar))^(ALP_PPI)*(exp(gdp)/gdp_ss)^ALP_GDP*exp(em);    //E28
// OptimalPolicy
//exp(ppiN)=BAR_PPI; 
// Crawling Peg
//exp(ppiS)=ppi_ss/ppist_ss ;

//FISCAL RULE AND TRANSFERS
exp(R(-1))*(d(-1))/exp(ppi)+exp(m)+exp(tao)=d+exp(m(-1))/exp(ppi(-1))+exp(pN)*exp(g);                     //E29
(exp(tao)-tao_ss)/gdp_ss = PHID*(((-1)*d+d_ss)/gdp_ss)+PHIG*(((exp(pN)*exp(g))-(exp(pN)*g_ss))/gdp_ss);      //E30

// Exogenous 
z-log(z_ss)            = RHO_z*(z(-1)-log(z_ss))+u_z;                     //E31
em-log(em_ss)          = RHO_em*(em(-1)-log(em_ss))+u_em;                 //E32
yT-log(yT_ss)          = RHO_yT*(yT(-1)-log(yT_ss))+u_yT;                 //E33
Rw-log(Rw_ss)          = RHO_Rw*(Rw(-1)-log(Rw_ss))+u_Rw;                 //E34
ppist-log(ppist_ss)    = RHO_ppist*(ppist(-1)-log(ppist_ss))+u_ppist;     //E35
ppiTar-log(BAR_PPI)    = RHO_ppiTar*(ppiTar(-1)-log(BAR_PPI))+u_ppiTar;   //E36
g-log(g_ss)            = RHO_g*(g(-1)-log(g_ss))+u_g;                     //E37
end; 

%----------------------------------------------------------------
% 4. Steady State
%----------------------------------------------------------------
steady_state_model; 


d_ss        = DEUDA;
ppi_ss      = BAR_PPI;
ppiTar_ss   = BAR_PPI;
ppiN_ss     = ppi_ss;
ppiT_ss     = ppi_ss;
Rw_ss       = 1/BET;
Rst_ss      = Rw_ss;
ppiS_ss     = ppi_ss/(Rst_ss*BET);
ppist_ss    = ppi_ss/ppiS_ss;

R_ss        = ppi_ss/BET;

mcN_ss      = (EPS-1)/(EPS)+PHIN/(EPS)*(ppiN_ss-1)*ppiN_ss*(1-BET);

yN_ss       = z_ss*h_ss^ALP;

g_ss        = sg_ss*(yT_ss+yN_ss);

cN_ss       = yN_ss*(1-PHIN/2*(ppiN_ss-1)^2)-g_ss;

dst_ss      = BARD;

dstC_ss     = dst_ss/(1-DELTA);  

tb_ss       = dst_ss*(Rst_ss/ppist_ss-1);

cT_ss       = yT_ss-tb_ss;

c_ss        = (GAM^(1/VRHO)*cN_ss^(1-1/VRHO) + (1-GAM)^(1/VRHO)*cT_ss^(1-1/VRHO))^(VRHO/(VRHO-1));

rer_ss      = (cT_ss/(c_ss*(1-GAM)))^(-1/VRHO);

p_ss        = 1/rer_ss;

pN_ss       = p_ss*(cN_ss/(GAM*c_ss))^(-1/VRHO);

w_ss        = pN_ss*rer_ss*mcN_ss*ALP*h_ss^(ALP-1);

hW_ss       = h_ss/DELTA; 

tao_ss      = d_ss*(1-R_ss)+pN_ss*g_ss;

cW_ss       =w_ss*hW_ss-tao_ss;

laW_ss      =cW_ss^(-SIG);

CHI         =laW_ss*w_ss/(hW_ss^VPHI);

cC_ss       =(c_ss-(DELTA*cW_ss))/(1-DELTA);

laC_ss      = cC_ss^(-SIG);

gdp_ss      = c_ss + yT_ss - cT_ss + g_ss;

stb_ss      = tb_ss*(c_ss/rer_ss+tb_ss);

syT_ss      = yT_ss*(c_ss/rer_ss+tb_ss);

mW_ss       = sm_ss*(c_ss+tb_ss*rer_ss);

m_ss        = DELTA*mW_ss;   

NU          = laW_ss*(1-BET/ppi_ss)*mW_ss^XI;

ngdp_ss     =(c_ss+tb_ss*rer_ss);

dC_ss       =d_ss/(1-DELTA);

//dC_ss       = d_ss/(1-DELTA);
//37 ecuaciones


//Definición para las 37 variables - hay 38 ecuaciones. 
c=log(c_ss);
laW=log(laW_ss);
laC=log(laC_ss);
h=log(h_ss);
w=log(w_ss);
Rst=log(Rst_ss);
ppiS=log(ppiS_ss);
ppi=log(ppi_ss);
R=log(R_ss);
cN=log(cN_ss);
pN=log(pN_ss);
cT=log(cT_ss);
ppiT=log(ppiT_ss);
rer=log(rer_ss);
ppiN=log(ppiN_ss);
mcN=log(mcN_ss);
yN=log(yN_ss);
gdp=log(gdp_ss);
dst=dst_ss;
tb=(tb_ss);
ppist=log(ppist_ss);
em=log(em_ss);
Rw=log(Rw_ss);
yT=log(yT_ss);
g=log(g_ss);      
ppiTar =log(ppiTar_ss);
m=log(m_ss);
z=log(z_ss);
p=log(p_ss);
hW=log(hW_ss);
cW=log(cW_ss);
cC=log(cC_ss);
mW=log(mW_ss);
tao=log(tao_ss);
d=d_ss;
dC=dC_ss;
dstC=dstC_ss;

end;
%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
//resid;
model_diagnostics;
steady;
check;

%----------------------------------------------------------------
% 7. Estimation
%----------------------------------------------------------------
shocks;
var u_z           = SIGMA_z^2;
var u_ppist       = SIGMA_ppist^2;
var u_em          = SIGMA_em^2;
var u_Rw          = SIGMA_Rw^2;
var u_yT          = SIGMA_yT^2;
var u_ppiTar      = SIGMA_ppiTar^2;
var u_g           = SIGMA_g^2;

end;

// RUNNING
stoch_simul(periods=0, irf=1000, order=1, nograph, ar=50);

save TNT_NK_taylor_punto2a.mat oo_ M_ options_;


