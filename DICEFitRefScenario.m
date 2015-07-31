function [Y,L,E,K,eland,forcoth,al,sigma,outyears]=DICEFitRefScenario(years,GDP,Population,IndustrialCO2,LandCO2,nonCO2forcings,tims)

% [Y,L,E,K,eland,forcoth,al,sigma]=DICEFitRefScenario(years,GDP,Population,IndustrialCO2,LandCO2,nonCO2forcings)
%
% Output DICE-appropriate parameters from an exogenously specified reference scenario.
%
% years: Years
% GDP: GDP in trillion 2005 dollars
% Population: Population in millions
% IndustrialCO2: Industrial and fossil CO2 emissions in Gt C/y
% LandCO2: Land use CO2 emissions in Gt C/y
% nonCO2forcings: non-CO2 forcing in W/m^2
%
% Example: MiniCAM-based scenarios for SCC analysis
%	
%	defp=DICEParameters;
%	years = 2000:10:2300;
%	GDP = [36.1 47.4 60.8 78.9 100.6 125.7 160.6 201.8 249.3 309.4 369.5 437.8 514.6 600.2 694.5 797.3 908.2 1026.3 1150.6 1279.8 1412.4 1543.4 1670.0 1789.3 1898.3 1994.2 2074.4 2136.7 2179.2 2200.8 2200.8];
%	%GDP2005 = (GDP(2)/GDP(1))^.5*GDP(1);
%	%GDP = GDP * defp.q0/GDP2005;
%	Population = [6.04 6.80 7.53 8.07 8.50 8.81 8.98 9.04 9.01 8.84 8.66 8.51 8.38 8.27 8.17 8.09 8.03 7.98 7.95 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93]*1000;
%	IndustrialCO2 = [26.53 31.81 38.01 45.14 51.74 57.80 62.75 67.78 72.89 76.67 80.46 83.77 86.53 88.68 90.17 90.96 91.03 90.39 89.05 87.04 84.40 81.04 77.06 72.55 67.63 62.43 57.06 51.65 46.28 41.07 36.09] * 12/44;
%	LandCO2 = [3.82 2.80 1.93 1.84 1.77 1.73 1.62 1.32 0.85 0.73 0.62 0.56 0.50 0.43 0.37 0.31 0.25 0.19 0.12 0.06 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
%	nonCO2forcings = [0.67 0.75 0.85 0.94 1.03 1.13 1.19 1.25 1.31 1.34 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41];
%	
%	[Y,L,E,K,eland,forcoth,al,sigma]=DICEFitRefScenario(years,GDP,Population,IndustrialCO2,LandCO2,nonCO2forcings);
%	
%	spec={'L',L,'al',al,'sigma',sigma,'etree',eland,'forcoth',forcoth};
%	p=DICEParameters('prstp',0.03,'elasmu',0,spec{:});
%	Ref=KDICE(spec{:},'aa2',0);
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Mar 4 18:20:03 EST 2013
% Copyright (C) 2013 by Robert E. Kopp; distributed under GNU GPL v3

defp = DICEParameters;
defval('tims',defp.t);
if max(years)>max(tims)
	deltat = tims(2)-tims(1);
	tims = [tims (tims(end)+deltat):deltat:(max(years)+deltat)];
end
outyears=tims;


Y = exp(interp1(years,log(GDP+1e-9),tims,'linear','extrap'));

L = exp(interp1(years,log(Population+1e-9),tims,'linear','extrap'));

E = exp(interp1(years,log(IndustrialCO2+1e-9),tims,'linear','extrap'));

eland = exp(interp1(years,log(LandCO2+1e-9),tims,'linear','extrap'))*10;

forcoth = exp(interp1(years,log(nonCO2forcings+1e-9),tims,'linear','extrap')); 


%--------------------------------------------------------------------------

% Solve for implied path of exogenous technical change

al     = zeros(1,length(tims));

al(1)  = defp.a0;

savings = 0.2;
K(1) = (Y(1)/al(1)/(L(1)^(1-defp.gama)))^(1/defp.gama);

for t = 2:length(tims);

    K(t) = K(t-1)*(1-defp.dk)^(tims(t)-tims(t-1)) + savings*Y(t-1)*(tims(t)-tims(t-1));

    al(t) = Y(t)/(K(t)^defp.gama * L(t)^(1-defp.gama));
end;
sigma = E./Y/(1-defp.miu_2005);


%--------------------------------------------------------------------------

