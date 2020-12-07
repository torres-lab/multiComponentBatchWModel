function [pHcalc,aHcalc,AlfreeCalc] = pHfromModel(Na,Ca,Al,pCO2)
% Aluminmum Dissociation Constants
k1 =10^(-4.987);
k2 =10^(-10.13);
k3 =10^(-16.76);
k4 =10^(-22.16);
% Carbon Dissociation Constants
kH =10^(-1.47);
ki =10^(-6.35);
aHinv = (Na + (2.*Ca) ) ./ (ki .* kH .* pCO2);
aH = 1./aHinv;
AlFree = Al./(1 + (k1./aH) + (k2./aH.^2) + (k3./aH.^3) + (k4./aH.^4));
AlNetCharge = AlFree.*(3 + (2.*k1./aH) + (k2./aH.^2) - (k4./aH.^4));
aHinv = (Na + (2.*Ca) +AlNetCharge) ./ (ki .* kH .* pCO2);
aHcalc = 1./aHinv;
pHcalc = -log10(aHcalc);
AlfreeCalc = Al./(1 + (k1./aHcalc) + (k2./aHcalc.^2) + (k3./aHcalc.^3) + (k4./aHcalc.^4));
