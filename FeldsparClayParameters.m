%% Feldspar Dissolution
ConversionFactor = 1E4; %cm2/m2 (rate law is norm to cm2)
fsparM = 1; %grams
fsparSA = 2; %reactive surface area (m2/gram)
temk = 3; %stoichiometric coefficient
mExp = 1; %near-Eq rate dependance (0.1)
nExp = 1; %near-Eq rate dependance (2)
DGo_albite = -11.2; 
DGo_anorthite = -137.7;
albRLcoef = [0.047, -0.6, 13.82; 0.049, -0.604, 13.25; 0.047, -0.6, 12.68];
anRLcoef = [0.09, -1.31, 10.4; 0.091, -1.36, 9.66; 0.085, -1.31, 9.2];
% Rate Law Coefficients
al1 = albRLcoef(rateLawModel,1);
al2 = albRLcoef(rateLawModel,2);
al3 = albRLcoef(rateLawModel,3);
an1 = anRLcoef(rateLawModel,1);
an2 = anRLcoef(rateLawModel,2);
an3 = anRLcoef(rateLawModel,3);
%% Clay Precipitation
kaolSA = 10; %reactive surface area (m2/gram)
KP = 2.21E-13; % RateConstant (mols/m2/sec) - Yang&Steefel2008
%% Aqueous Speciation
% Aluminmum Dissociation Constants
k1 =10^(-4.987);
k2 =10^(-10.13);
k3 =10^(-16.76);
k4 =10^(-22.16);
% Carbon Dissociation Constants
kH =10^(-1.47);
ki =10^(-6.35);
%% Model Boundary / Initial Conditions + Constants
Rcnst = 8.314./1000; %kJ K-1 mol-1
Temp = 298; %kelvin (dont change as model not built for multiple temps)
options = odeset('reltol',1e-5,'abstol',1e-9);