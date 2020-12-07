%% Set Rate Law
rateLawModel = 2; %slow (1), best-fit (2), or fast (3) rate law
run('FeldsparClayParameters.m')
PrcntAn = 0.5; % percent of anorthite in feldspar
PrcntAl = 1-PrcntAn; % percent of albite in feldspar
%% Set pCO2
pCO2 = 0.01; %bars
%% Set Clay Species
DGo_kaolinite = -42.88; %-42.88 = HALLOYSITE | -23.63 = KAOLINITE
kaolM = 0.1; % initial mass of Kaolinite (g)
%% set forced Al
Alslope = 4E-14; %Al concentration slope for linear Al model
Alval = 0.5E-6; %Al concentration for constant Al model
%% Set model initial conditions and duration
x0 = [1E-6,1E-6,1E-6,1E-12]; %initial concentrations (Na, Ca, Si, Al; molar)
tlengthYears = 1.5; %simulation length in years
tlength = tlengthYears*365*24*60*60; %simulation length in seconds
%% Set W/R
WR = 0.9; %water/rock ratio
fVol = (fsparM.*fsparSA).*WR; % fluid volume from W/R (liters)
%% ODE solver
modelOption = 4;
if modelOption == 1 %variable [Al], linear precip rate law
[T,C] = ode23t(@varAlLinP,[1,tlength],x0(1:4),options,fVol,Rcnst,Temp,pCO2,...
    fsparM,fsparSA,temk,mExp,nExp,DGo_albite,DGo_anorthite,PrcntAn,PrcntAl,...
    kaolM,kaolSA,KP,DGo_kaolinite,...
    k1,k2,k3,k4,kH,ki,ConversionFactor,an1,an2,an3,al1,al2,al3);
    % Calculate Al3+ and pH from model output 
    [pHc,aHc,Al3c] = pHfromModel(C(:,1),C(:,2),C(:,4),pCO2);

elseif modelOption == 2 %linear [Al], linear precip rate law
[T,C] = ode23t(@linAlLinP,[1,tlength],x0(1:3),options,fVol,Rcnst,Temp,pCO2,...
    fsparM,fsparSA,temk,mExp,nExp,DGo_albite,DGo_anorthite,PrcntAn,PrcntAl,...
    kaolM,kaolSA,KP,DGo_kaolinite,...
    k1,k2,k3,k4,kH,ki,ConversionFactor,an1,an2,an3,al1,al2,al3,Alslope);
    % Calculate Al3+ and pH from model output 
    [pHc,aHc,Al3c] = pHfromModel(C(:,1),C(:,2),((Alslope).*T)+1E-12,pCO2);

elseif modelOption == 3 %constant [Al], linear precip rate law
[T,C] = ode23t(@conAlLinP,[1,tlength],x0(1:3),options,fVol,Rcnst,Temp,pCO2,...
    fsparM,fsparSA,temk,mExp,nExp,DGo_albite,DGo_anorthite,PrcntAn,PrcntAl,...
    kaolM,kaolSA,KP,DGo_kaolinite,...
    k1,k2,k3,k4,kH,ki,ConversionFactor,an1,an2,an3,al1,al2,al3,Alval);
    % Calculate Al3+ and pH from model output 
    [pHc,aHc,Al3c] = pHfromModel(C(:,1),C(:,2),Alval,pCO2);
    
elseif modelOption == 4 %variable [Al], TST precip rate law
[T,C] = ode23t(@varAlTSTP,[1,tlength],x0(1:4),options,fVol,Rcnst,Temp,pCO2,...
    fsparM,fsparSA,temk,mExp,nExp,DGo_albite,DGo_anorthite,PrcntAn,PrcntAl,...
    kaolM,kaolSA,KP,DGo_kaolinite,...
    k1,k2,k3,k4,kH,ki,ConversionFactor,an1,an2,an3,al1,al2,al3);
    % Calculate Al3+ and pH from model output 
    [pHc,aHc,Al3c] = pHfromModel(C(:,1),C(:,2),C(:,4),pCO2);
    
elseif modelOption == 5 %constant [Al], TST precip rate law
[T,C] = ode23t(@conAlTSTP,[1,tlength],x0(1:3),options,fVol,Rcnst,Temp,pCO2,...
    fsparM,fsparSA,temk,mExp,nExp,DGo_albite,DGo_anorthite,PrcntAn,PrcntAl,...
    kaolM,kaolSA,KP,DGo_kaolinite,...
    k1,k2,k3,k4,kH,ki,ConversionFactor,an1,an2,an3,al1,al2,al3,Alval);
    % Calculate Al3+ and pH from model output 
    [pHc,aHc,Al3c] = pHfromModel(C(:,1),C(:,2),Alval,pCO2);
    
else
end

%% Calculate saturation from model output 
Qkao = ((Al3c.^2).*(C(:,3).^2))./(aHc.^6); %reaction quotient
deltaGkao = DGo_kaolinite + (Rcnst.*Temp.*log(Qkao)); %delta G
dGtermKao = deltaGkao./(1.*Rcnst.*Temp); %deltaG/RT

%% Example Plot
figure
subplot(1,3,1)
hold on
plot(T./60./60./24./365,(C(:,1)+C(:,2)).*1E6,'-','linewidth',2)
xlabel('Time (years)'); ylabel('Na+Ca (\muM)')
subplot(1,3,2)
hold on
plot(T./60./60./24./365,C(:,3).*1E6,'-','linewidth',2)
xlabel('Time (years)'); ylabel('Si (\muM)')
subplot(1,3,3)
hold on
plot(T./60./60./24./365,pHc,'-','linewidth',2)
xlabel('Time (years)'); ylabel('pH')
