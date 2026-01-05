function dydt = EDC_ode(t, y, option, param)


%% -------------------------- PARAMETERS MAPPING ----------------------------------%%

%Manual code but fast to execute

k1 = param.k1;
k2 = param.k2;
k3 = param.k3;
k30 = param.k30;
wc = param.wc;
Kd3 = param.Kd3;
n3 = param.n3;
k4 = param.k4;

kf7 = param.kf7;
kb7 = param.kb7;
kf8 = param.kf8; 
kb8 = param.kb8; 
kf5 = param.kf5;
kb5 = param.kb5;
kf6 = param.kf6;
kb6 = param.kb6;
X = param.X;
CRtot = param.CRtot; 
PRtot = param.PRtot; 

%% ------------------------- STATE NAME MAPPING----------------------------%%

EH   = y(1);
PH   = y(2);
EHCR = y(3);
XCR  = y(4);
EHPR = y(5);
XPR  = y(6);

CR = CRtot - EHCR - XCR;
PR = PRtot - EHPR - XPR;

%% ------------------------------ ODEs-------------------------------------%%

dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%EH
dydt(1) = k1*PH - k2*EH;

%PH
dydt(2) = k30 + k3*Kd3^n3/(Kd3^n3+(EHCR+wc*XCR)^n3) - k4*PH;

%EHCR
dydt(3) = kf7*EH*CR - kb7*EHCR;

%XCR
dydt(4) = kf8*X*CR - kb8*XCR;

%EHPR
dydt(5) = kf5*EH*PR - kb5*EHPR;

%XPR
dydt(6) = kf6*X*PR - kb6*XPR;

end