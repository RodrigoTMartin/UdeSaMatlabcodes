%% TNT SGU - Rotemeberg prices

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

// Endogenous variables 
var c la h w Rst ppiS ppi R cN pN p cT ppiT rer ppiN mcN yN m ;
var gdp dst tb;

// Variables that hav associated an exogenous process
var z ppist em Rw yT ppiTar;

// Disturbances variables (7)
varexo u_z u_ppist u_em u_Rw u_yT u_ppiTar;

// Parameters 
parameters SIG VPHI CHI BET GAM VRHO NU XI;
parameters ALP EPS THET PHIN;
parameters ALP_PPI ALP_GDP BAR_PPI PHI_B BARD;
parameters RHO_z RHO_ppist RHO_em RHO_Rw RHO_yT RHO_ppiTar;
parameters SIGMA_z SIGMA_ppist SIGMA_em SIGMA_Rw SIGMA_yT SIGMA_ppiTar;

parameters R_ss gdp_ss ppi_ss z_ss ppist_ss em_ss Rw_ss yT_ss ngdp_ss; // in eqns and exon
parameters h_ss stb_ss rer_ss sm_ss syT_ss;                      // given in ss


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
// Calibarted parameteres and steady state values
% CALCULATED ENDOGENOUSLY
% ppist_ss, CHI, BARD, delF_ss, yT_ss

BET          = 1/1.03;   // SGU
SIG          = 2;             // SGU
VPHI         = 0.5;           // SGU
VRHO         = 0.5;           // SGU
GAM          = 1-0.26;        // SGU

ALP        = 1;          // SGU     
EPS        = 6;             // SGU
THET       = 0.7;           // SGU

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

SIGMA_z      = 0.01;
SIGMA_ppist  = 0.01;
SIGMA_em     = 0.01;
SIGMA_Rw     = 0.01;
SIGMA_yT     = 0.01;
SIGMA_ppiTar = 0.01;

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




%----------------------------------------------------------------
% 3. Model 
%----------------------------------------------------------------
model; 
// COMMON FOR ALL MODELS
// HH 7 eqns
exp(la) = exp(c)^(-SIG);                                             // E1
CHI*exp(h)^VPHI = exp(la)*exp(w);                                    // E2 
exp(la) = BET*exp(R)*exp(la(+1))/exp(ppi(+1));                       // E4
exp(la) = BET*exp(Rst)*exp(la(+1))*exp(ppiS(+1))/exp(ppi(+1));       // E3
exp(la) = BET*exp(la(+1))/exp(ppi(+1)) + NU*exp(m)^(-XI);                       // E4

exp(cN) = GAM*(exp(pN)/exp(p))^(-VRHO)*exp(c);                                // E5
exp(cT) = (1-GAM)*(1/exp(p))^(-VRHO)*exp(c);                                                                          // E6
exp(c) = (GAM^(1/VRHO)*exp(cN)^(1-1/VRHO) + (1-GAM)^(1/VRHO)*exp(cT)^(1-1/VRHO))^(VRHO/(VRHO-1));


exp(yN)=exp(z)*exp(h)^(ALP);

exp(pN)*exp(rer)*exp(mcN)*exp(z)*ALP*exp(h)^(ALP-1)=exp(w);


EPS-1=exp(mcN)*(EPS)-PHIN*(exp(ppiN)-1)*exp(ppiN)+PHIN*BET*exp(la(+1))/(exp(ppi(+1))*exp(la))*exp(yN(+1))/exp(yN)*(exp(ppiN(+1))-1)*(exp(ppiN(+1)))^2;


// Foreign Sector 7 eqns
exp(rer)/exp(rer(-1))= exp(ppiS)*exp(ppist)/exp(ppi);                                                           // E14

exp(rer)=1/exp(p);
exp(Rst) = exp(Rw)*exp(PHI_B*(dst-BARD));                                                           // E15

// Equilibrium conditions and closing model 1 eqns
exp(cN) = exp(yN)*(1-PHIN/2*(exp(ppiN)-1)^2);

// National identities and definitions 4 eqns
exp(Rst(-1))/exp(ppist)*dst(-1) = tb + dst;                                                           // E19

exp(gdp) = exp(c) + tb;                                                                                           // E21

exp(ppiT) = exp(rer)*exp(ppi)/exp(rer(-1));                                                                         // E11

exp(ppiN) = exp(pN)*exp(ppiT)/exp(pN(-1));   

tb = exp(yT)-exp(cT);                                                                                    // E18


//Monetary Policy
// Taylor Rule
exp(R)= R_ss*exp(ppiTar)*(exp(ppi)/exp(ppiTar))^(ALP_PPI)*(exp(gdp)/gdp_ss)^ALP_GDP*exp(em);            //E13

// OptimalPolicy
//exp(ppiN)=BAR_PPI; 

// Crawling Peg
//exp(ppiS)=ppi_ss/ppist_ss ;

// Exogenous 
z-log(z_ss)            = RHO_z*(z(-1)-log(z_ss))+u_z;             //E21
em-log(em_ss)          = RHO_em*(em(-1)-log(em_ss))+u_em;             //E21
yT-log(yT_ss)          = RHO_yT*(yT(-1)-log(yT_ss))+u_yT;             //E22
Rw-log(Rw_ss)          = RHO_Rw*(Rw(-1)-log(Rw_ss))+u_Rw;             //E23
ppist-log(ppist_ss)    = RHO_ppist*(ppist(-1)-log(ppist_ss))+u_ppist; //E24
ppiTar-log(BAR_PPI)    = RHO_ppiTar*(ppiTar(-1)-log(BAR_PPI))+u_ppiTar; //E24

end; 

%----------------------------------------------------------------
% 4. Steady State
%----------------------------------------------------------------
steady_state_model; 

ppi_ss      = BAR_PPI;
ppiTar_ss   = BAR_PPI;
ppiN_ss     = ppi_ss;
ppiT_ss     = ppi_ss;
Rw_ss       = 1/BET;
Rst_ss      = Rw_ss;
ppiS_ss     = ´ñppi_ss/(Rst_ss*BET);
ppist_ss    = ppi_ss/ppiS_ss;

mcN_ss      = (EPS-1)/(EPS)+PHIN/(EPS)*(ppiN_ss-1)*ppiN_ss*(1-BET);

yN_ss       = z_ss*h_ss^ALP;

cN_ss       = yN_ss*(1-PHIN/2*(ppiN_ss-1)^2);

dst_ss      = BARD;

tb_ss       = dst_ss*(Rst_ss/ppist_ss-1);

cT_ss       = yT_ss-tb_ss;

c_ss        = (GAM^(1/VRHO)*cN_ss^(1-1/VRHO) + (1-GAM)^(1/VRHO)*cT_ss^(1-1/VRHO))^(VRHO/(VRHO-1));

rer_ss      = (cT_ss/(c_ss*(1-GAM)))^(-1/VRHO);

p_ss        = 1/rer_ss;

pN_ss       = p_ss*(cN_ss/(GAM*c_ss))^(-1/VRHO);

la_ss       = c_ss^(-SIG);

w_ss        = pN_ss*rer_ss*mcN_ss*ALP*h_ss^(ALP-1);

CHI         = la_ss*w_ss/(h_ss^VPHI);

gdp_ss      = c_ss + yT_ss - cT_ss;

stb_ss      = tb_ss*(c_ss/rer_ss+tb_ss);

syT_ss      = yT_ss*(c_ss/rer_ss+tb_ss);

m_ss        = sm_ss*(c_ss+tb_ss*rer_ss);

NU          = la_ss*(1-BET/ppi_ss)*m_ss^XI;

ngdp_ss=(c_ss+tb_ss*rer_ss);

R_ss        = ppi_ss/BET;

c=log(c_ss);
la=log(la_ss);
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
dst=(dst_ss);
tb=(tb_ss);
ppist=log(ppist_ss);
em=log(em_ss);
Rw=log(Rw_ss);
yT=log(yT_ss);
ppiTar =log(ppiTar_ss);
m=log(m_ss);
z=log(z_ss);
p=log(p_ss);

end;
%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
steady;
check;

%----------------------------------------------------------------
% 7. Estimation
%----------------------------------------------------------------
shocks;
var u_z          = SIGMA_z^2;
var u_ppist       = SIGMA_ppist^2;
var u_em          = SIGMA_em^2;
var u_Rw          = SIGMA_Rw^2;
var u_yT          = SIGMA_yT^2;
var u_ppiTar       = SIGMA_ppiTar^2;

end;

// RUNNING
stoch_simul(periods=0, irf=1000, order=1, nograph, ar=50);

save TNT_NK_Taylor.mat oo_ M_ options_;


