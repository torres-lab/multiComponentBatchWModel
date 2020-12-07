function xprime = conAlTSTP(T,C,fVol,Rcnst,Temp,pCO2,...
    fsparM,fsparSA,temk,mExp,nExp,DGo_albite,DGo_anorthite,PrcntAn,PrcntAl,...
    kaolM,kaolSA,KP,DGo_kaolinite,...
    k1,k2,k3,k4,kH,ki,ConversionFactor,an1,an2,an3,al1,al2,al3,Alval)
%% Calculate pH from Na and Ca
aHinv = (C(1) + (2.*C(2)) ) ./ (ki .* kH .* pCO2);
aH = 1./aHinv;

%% Calculate Al3+ from pH and total Al
AlFree = Alval./...
    (1 + (k1./aH) + (k2./aH.^2) + (k3./aH.^3) + (k4./aH.^4));
AlNetCharge = AlFree.*(3 + (2.*k1./aH) + (k2./aH.^2) - (k4./aH.^4));
aHinv = (C(1) + (2.*C(2)) +AlNetCharge) ./ (ki .* kH .* pCO2);
aH = 1./aHinv;
pH = -log10(aH);
AlFree = Alval./...
    (1 + (k1./aH) + (k2./aH.^2) + (k3./aH.^3) + (k4./aH.^4));
%% Calculate Kaolinite Thermodynamics
Qkao = ((AlFree.^2).*(C(3).^2))./(aH.^6); %reaction quotient
deltaGkao = DGo_kaolinite + (Rcnst.*Temp.*log(Qkao)); %delta G
dGtermKao = deltaGkao./(2.*Rcnst.*Temp); %deltaG/RT

%% Calculate Plagioclase Thermodynamics 
Qalb = (C(1).*AlFree.*(C(3).^3))./(aH.^4); 
DGalb= DGo_albite + (Rcnst.*Temp.*log(Qalb));
DGTerm_alb = DGalb./(temk.*Rcnst.*Temp);
%
Qan = (C(2).*(AlFree.^2).*(C(3).^2))./(aH.^8);  
DGan= DGo_anorthite + (Rcnst.*Temp.*log(Qan));
DGTerm_an = DGan./(temk.*Rcnst.*Temp);

%% Dissolution Rate Laws
AnRateSiLog = @(pH) (an1 .* (pH.^2)) + (an2 .*pH) - an3; 
AnDissRate = (10.^(AnRateSiLog(pH))) .* ((1-(exp(DGTerm_an).^mExp)).^nExp);
AnDissFlux = PrcntAn.*fsparM.*fsparSA.*AnDissRate.*ConversionFactor.*heaviside(-DGan);
%
AlbRateSiLog = @(pH) (al1 .* (pH.^2)) + (al2.*pH) - al3; 
AlbDissRate = (10.^(AlbRateSiLog(pH))) .* ((1-(exp(DGTerm_alb).^mExp)).^nExp);
AlbDissFlux = PrcntAl.*fsparM.*fsparSA.*AlbDissRate.*ConversionFactor;

%% Precipitation Rate Law
precipRate = (KP./10) .* (exp(dGtermKao)-1);
precipFlux = heaviside(deltaGkao).*(2.*precipRate.*kaolM.*kaolSA);
%% Differential Equations
xprime = zeros(3,1); %Na, Ca, Si, Al
xprime(1) = (AlbDissFlux./3) ./ fVol;                                   %Na
xprime(2) = (AnDissFlux./2) ./ fVol;                                    %Ca
xprime(3) = ( (AlbDissFlux + AnDissFlux) - (precipFlux) ) ./ fVol;      %Si