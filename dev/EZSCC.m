function SCCout=EZSCC(p,SavingsRate,miu,timeslot,eta,RRA,rho)

% SCCout = EZSCC(p,[SavingsRate],[miu],[timeslot],[eta],[RRA],[rho])
%
% INPUTS
%
% p takes the form of DICEParameters
% SavingsRate is a 1 x t vector of savings rates (defaults to p.basesavings)
% miu is a 1 x t vector of abatement rates (defaults to p.miu_2005)
% timeslot is the index of the time slot in which the incremental changes will be made (defaults to 2)
% damagemodfactor is a multiplier of damages (defaults to 1)
%
% Values in SCCout are denominated per ton CO2.
%
% Example:
%
%   defp = DICEParameters;
%
%   cses = randRoeBaker([1000 1]);
%   aa3s = randTri([1000 1],[1 1.5 3.5]);
%   aa2s = defp.aa2 * 2.5.^(defp.aa3-aa3s);
%
%   spec={'elasmu',0,'prstp',.03,'T2xCO2',cses,'aa3',aa3s,'aa2',aa2s};
%   p=DICEParameters(spec{:});
%
%   pexp=DICEParameters(spec{:},'calcdamages','1 - exp(-p.aa2.*abs(Temperature).^p.aa3)')
%
%   SCCs = SCC(p);
%   SCCsexp = SCC(pexp);
%	
%	clf;
%	subplot(2,1,1)
%	hist([SCCs.SCC_consumptiondenominated],20);
%	xlabel('SCC ($/ton CO_2)');
%	ylabel('N');
%	title('Damages = 1 - 1/(1 + aT^b)')
%
%   subplot(2,1,2)
%	hist([SCCsexp.SCC_consumptiondenominated],20);
%	xlabel('SCC ($/ton CO_2)');
%	ylabel('N');
%	title('Damages = 1 - exp(-aT^b)')
%
% Last updated by Bob Kopp rkopp-at-alumni.caltech.edu, 20 May 2011

	defval('L',p.L);
	defval('t',p.t);
	defval('eta',1./p.elasmu);
	defval('RRA',p.elasmu);
	defval('rho',p.prstp);
	defval('damagemodfactor',1);
	defval('SavingsRate',p.basesavings);
	defval('miu',ones(size(p.t))*p.miu_2005);
	defval('timeslot',2);
	
	SCCout.eta = eta;
	SCCout.RRA = RRA;
	SCCout.rho = rho;
	
	
	EZU = @(C) EZUtility(C,L,eta,RRA,rho,t);
	
	N = length(p.T2xCO2);
	
	t=p.t;
	incrementalemissions = zeros(size(t));
	incrementalconsumption = zeros(size(t));
	incrementalcapital = zeros(size(t));
	
	incrementalemissions(timeslot) = 1e-3; % = 1 million tonnes
	incrementalconsumption(timeslot) = 1e-3; % = 1 billion dollars
	incrementalcapital(timeslot) = 1e-3; % = 1 billion dollars

	%convfactor = @(scen,p,ecosh,ecoel) ( (scen.ConsumptionPerCapita./repmat(scen.ConsumptionPerCapita(:,1),1,length(p.t))).^(1./ecoel) .*  repmat(scen.EffectiveConsumptionPerCapita(:,1),1,length(p.t))./scen.EffectiveConsumptionPerCapita  ) .^ (1 - p.elasmu .* ecoel) ;

	[EW,baseline] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor);		
	[EW,plusemissions] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor,incrementalemissions);
	[EW,plusconsumption] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor,zeros(size(t)),incrementalconsumption);
	[EW,pluscapital] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor,zeros(size(t)),zeros(size(t)),incrementalcapital);
	
	ecoel=p.ecoelasticity; ecosh = p.ecoshare;
	if length(ecoel)>1
		ecoel = repmat(ecoel,1,length(p.t));
	end
	if length(ecosh)>1
		ecosh = repmat(ecosh,1,length(p.t));
	end
	%convfactors = real(convfactor(baseline,p,ecosh,ecoel));
	%convfactors = (baseline.ConsumptionPerCapita./baseline.EffectiveConsumptionPerCapita)./repmat((baseline.ConsumptionPerCapita(:,1)./baseline.EffectiveConsumptionPerCapita(:,1)),1,length(p.t));
	
	
	deltaemissions = (EZU(plusemissions.EffectiveConsumptionPerCapita)-EZU(baseline.EffectiveConsumptionPerCapita))/sum(incrementalemissions);
	deltaconsumption = (EZU(plusconsumption.EffectiveConsumptionPerCapita)-EZU(baseline.EffectiveConsumptionPerCapita))/sum(incrementalconsumption);
	deltacapital = (EZU(pluscapital.EffectiveConsumptionPerCapita)-EZU(baseline.EffectiveConsumptionPerCapita))/sum(incrementalcapital);
	SCCout.ESCC_consumptiondenominated = (12/44)*-deltaemissions./deltaconsumption*1e3;
	%SCCout.ESCC_investmentdenominated = (12/44)*-deltaemissions./deltacapital*1e3;

%	SCCout.deltaemissions = deltaemissions;
%	SCCout.deltaconsumption = deltaconsumption;
%	SCCout.deltacapital = deltacapital;
%
	%convemissions = sum((plusemissions.UtilityPerCapita - baseline.UtilityPerCapita) .* convfactors .* repmat(p.L .* p.rr,N,1),2);
	%convconsumption = sum((plusconsumption.UtilityPerCapita - baseline.UtilityPerCapita) .* convfactors .* repmat(p.L .* p.rr,N,1),2);

%	convemissions = (plusemissions.EquivalentMaterialWelfare-baseline.EquivalentMaterialWelfare)/sum(incrementalemissions);
%	convconsumption = (plusconsumption.EquivalentMaterialWelfare-baseline.EquivalentMaterialWelfare)/sum(incrementalemissions);
%
%	SCCout.SCC_materialconsumptionequiv = convemissions./convconsumption * (-12/44) * sum(incrementalcapital)/sum(incrementalemissions)*1e3;
%
%	incrementalabatement = zeros(size(t));
%	incrementalabatement(timeslot) = 1e-3; % = 1 million tonnes
%	incrementalmiu = repmat(incrementalabatement,N,1) ./ ((p.t(2)-p.t(1)).*repmat(p.sigma,N,1).*baseline.Output_Gross+repmat(p.etree,N,1));
%	incrementalforcing = repmat(p.forcoth,N,1) .* incrementalmiu;
%	
%	[EW,plusabatement] = DICEEconomicModel(p,SavingsRate,repmat(miu,N,1)+incrementalmiu,damagemodfactor,incrementalabatement,[],[],incrementalforcing);
%
%	deltaabatement = (plusabatement.Welfare-baseline.Welfare)/sum(incrementalabatement);
%	SCCout.MAC_consumptiondenominated = (12/44)*-deltaabatement./deltaconsumption*1e3;
%	SCCout.MAC_investmentdenominated = (12/44)*-deltaabatement./deltacapital*1e3;
%
%	Edeltaemissions = (plusemissions.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalemissions);
%	Edeltaconsumption = (plusconsumption.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalconsumption);
%	Edeltacapital = (pluscapital.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalcapital);
%	Econvemissions = (plusemissions.EquivalentMaterialExpectedWelfare-baseline.EquivalentMaterialExpectedWelfare)/sum(incrementalemissions);
%	Econvconsumption = (plusconsumption.EquivalentMaterialExpectedWelfare-baseline.EquivalentMaterialExpectedWelfare)/sum(incrementalemissions);
%	SCCout.ESCC_consumptiondenominated = (12/44)*-Edeltaemissions./Edeltaconsumption*1e3;
%	SCCout.ESCC_investmentdenominated = (12/44)*-Edeltaemissions./Edeltacapital*1e3;
%	SCCout.ESCC_materialconsumptionequiv = (12/44)*-Econvemissions./Econvconsumption*1e3;;
%
%	Edeltaabatement = (plusabatement.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalabatement);
%	SCCout.EMAC_consumptiondenominated = (12/44)*-Edeltaabatement./Edeltaconsumption*1e3;
%	SCCout.EMAC_investmentdenominated = (12/44)*-Edeltaabatement./Edeltacapital*1e3;

end
