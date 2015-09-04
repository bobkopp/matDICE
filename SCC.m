function [SCCout,scenarios]=SCC(p,SavingsRate,miu,timeslot,damagemodfactor,altmiu,liability)

% SCCout = SCC(p,[SavingsRate],[miu],[timeslot],[damagemodfactor],[altmiu],[liability])
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
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jul 26 23:29:21 EDT 2012
% Copyright (C) 2012 by Robert E. Kopp; distributed under GNU GPL v3

	defval('damagemodfactor',1);

	defval('SavingsRate',p.basesavings);
	defval('miu',ones(size(p.t))*p.miu0);
	defval('timeslot',2);
	defval('altmiu',miu);
	defval('liability',0);
	
	if p.doROI == 0
		liability = 0;
	end
	
	N = length(p.T2xCO2);
	

	t=p.t;


	incrementalemissions = zeros(size(t));

	incrementalconsumption = zeros(size(t));

	incrementalcapital = zeros(size(t));

	

	incrementalemissions(timeslot) = 1e-2; % = 10 million tonnes

	incrementalconsumption(timeslot) = 1e-2; % = 10 billion dollars

	incrementalcapital(timeslot) = 1e-2; % = 10 billion dollars


	%convfactor = @(scen,p,ecosh,ecoel) ( (scen.ConsumptionPerCapita./repmat(scen.ConsumptionPerCapita(:,1),1,length(p.t))).^(1./ecoel) .*  repmat(scen.EffectiveConsumptionPerCapita(:,1),1,length(p.t))./scen.EffectiveConsumptionPerCapita  ) .^ (1 - p.elasmu .* ecoel) ;

	[EW,ref] = DICEEconomicModel(p,SavingsRate,miu*0,damagemodfactor);		

	[EW,baseline] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor);		

	[EW,plusemissions] = DICEEconomicModel(p,SavingsRate,altmiu,damagemodfactor,incrementalemissions);
	
	if liability ~= 0
		[EW,plusemissions] = DICEEconomicModel(p,SavingsRate,altmiu,damagemodfactor,incrementalemissions,liability*(t(2)-t(1))*((plusemissions.Consumption-baseline.Consumption)-(plusemissions.ROI.Consumption-baseline.ROI.Consumption)));
	end	
	
	[EW,plusconsumption] = DICEEconomicModel(p,SavingsRate,altmiu,damagemodfactor,zeros(size(t)),incrementalconsumption);

	[EW,pluscapital] = DICEEconomicModel(p,SavingsRate,altmiu,damagemodfactor,zeros(size(t)),zeros(size(t)),incrementalcapital);

	
	ecoel=p.ecoelasticity; ecosh = p.ecoshare;
	if length(ecoel)>1
		ecoel = repmat(ecoel,1,length(p.t));
	end
	if length(ecosh)>1
		ecosh = repmat(ecosh,1,length(p.t));
	end
	%convfactors = real(convfactor(baseline,p,ecosh,ecoel));
	%convfactors = (baseline.ConsumptionPerCapita./baseline.EffectiveConsumptionPerCapita)./repmat((baseline.ConsumptionPerCapita(:,1)./baseline.EffectiveConsumptionPerCapita(:,1)),1,length(p.t));
	

	deltaemissions = (plusemissions.Welfare-baseline.Welfare)/sum(incrementalemissions);

	deltaconsumption = (plusconsumption.Welfare-baseline.Welfare)/sum(incrementalconsumption);

	deltacapital = (pluscapital.Welfare-baseline.Welfare)/sum(incrementalcapital);

	SCCout.SCC_consumptiondenominated = (12/44)*-deltaemissions./deltaconsumption*1e3;

	SCCout.SCC_investmentdenominated = (12/44)*-deltaemissions./deltacapital*1e3;


	%convemissions = sum((plusemissions.UtilityPerCapita - baseline.UtilityPerCapita) .* convfactors .* repmat(p.L .* p.rr,N,1),2);
	%convconsumption = sum((plusconsumption.UtilityPerCapita - baseline.UtilityPerCapita) .* convfactors .* repmat(p.L .* p.rr,N,1),2);

	convemissions = (plusemissions.EquivalentMaterialWelfare-baseline.EquivalentMaterialWelfare)/sum(incrementalemissions);
	convconsumption = (plusconsumption.EquivalentMaterialWelfare-baseline.EquivalentMaterialWelfare)/sum(incrementalemissions);

	SCCout.SCC_materialconsumptionequiv = convemissions./convconsumption * (-12/44) * sum(incrementalcapital)/sum(incrementalemissions)*1e3;

	incrementalmiu = zeros(size(t));
	incrementalmiu(timeslot) = 1e-3;
	incrementalmiu = repmat(incrementalmiu,N,1);
		
	if size(miu,1)<N
		miu = repmat(miu,N,1);
	end
	
	[EW,plusabate0] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor,[],[],[],[],incrementalmiu);
	incrementalforcing = baseline.Forcing-plusabate0.Forcing;
	[EW,plusabatement] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor,[],[],[],incrementalforcing,incrementalmiu);

 calcmargabatecost = @(PartFraction,CostFactor,expcost2,Mitigation,output,emit0)  bsxfun(@times,PartFraction.^(1-expcost2).*expcost2.*CostFactor,Mitigation.^(expcost2-1) .*(1e3*output*10)./(emit0*(44/12)));


	if p.doROI
		incrementalabatement = (baseline.ROI.Emissions - plusabatement.ROI.Emissions); % this is in GtC/time unit
		deltaabatement = (plusabatement.Welfare-baseline.Welfare)./sum(incrementalabatement,2);
		SCCout.MAC_consumptiondenominated = (12/44)*-deltaabatement./deltaconsumption*1e3;
		SCCout.MAC_investmentdenominated = (12/44)*-deltaabatement./deltacapital*1e3;
		SCCout.MAC0 = bsxfun(@times,(12/44)*((plusabatement.ROI.AbatementCost(:,timeslot).*plusabatement.ROI.Output_Gross(:,timeslot))-(baseline.ROI.AbatementCost(:,timeslot).*baseline.ROI.Output_Gross(:,timeslot)))*(t(2)-t(1))./(incrementalabatement(:,timeslot)+eps) * 1e3,(1-SavingsRate));
		%SCCout.MAC1 = calcmargabatecost(p.partfract,p.cost1,p.expcost2,miu,baseline.Output_Gross,baseline.Emissions_Gross);
		SCCout.MAC1=bsxfun(@times,calcmargabatecost(1,p.cost1,p.expcost2,baseline.ROI.Abatement,baseline.ROI.Output_Gross,baseline.ROI.Emissions_Gross),(1-SavingsRate));

		SCCout.EMAC1 = mean(SCCout.MAC1,1);
	else
		incrementalabatement = (baseline.Emissions - plusabatement.Emissions); % this is in GtC/time unit
		deltaabatement = (plusabatement.Welfare-baseline.Welfare)./sum(incrementalabatement,2);
		SCCout.MAC_consumptiondenominated = (12/44)*-deltaabatement./deltaconsumption*1e3;
		SCCout.MAC_investmentdenominated = (12/44)*-deltaabatement./deltacapital*1e3;
		SCCout.MAC0 = bsxfun(@times,(12/44)*((plusabatement.AbatementCost(:,timeslot).*plusabatement.Output_Gross(:,timeslot))-(baseline.AbatementCost(:,timeslot).*baseline.Output_Gross(:,timeslot)))*(t(2)-t(1))./(incrementalabatement(:,timeslot)+eps) * 1e3, 1-SavingsRate);
		SCCout.MAC1 = bsxfun(@times,calcmargabatecost(p.partfract,p.cost1,p.expcost2,miu,baseline.Output_Gross,baseline.Emissions_Gross),1-SavingsRate);
		SCCout.EMAC1 = mean(SCCout.MAC1,1);
	end

	baseline.AbatedCEmissions = baseline.abatement.*bsxfun(@plus,(t(2)-t(1))*bsxfun(@times,p.sigma,baseline.Output_Gross),p.etree);
	baseline.NonCO2ForcingReduction = bsxfun(@times,baseline.abatementnonco2,p.forcoth);
	AbatedCEmissionsPerMiu = baseline.AbatedCEmissions + baseline.Emissions;
	nonco2miupermiu = ( (baseline.abatementnonco2<p.forcothlimmiu)) * p.forcothmiufactor;

	effdeltamiu = bsxfun(@times,incrementalemissions,1./(eps+AbatedCEmissionsPerMiu));
	effdeltanonco2miu = zeros(size(effdeltamiu));
	effdeltanonco2miu(:,2:end) = [effdeltamiu(:,1:end-1).*nonco2miupermiu(:,2:end)];
	incrementalnonCO2 = bsxfun(@times,effdeltanonco2miu,p.forcoth);
	[EW,plusemissionsplusnonCO2] = DICEEconomicModel(p,SavingsRate,altmiu,damagemodfactor,incrementalemissions,[],[],incrementalnonCO2);
	
	if liability ~= 0
		[EW,plusemissionsplusnonCO2] = DICEEconomicModel(p,SavingsRate,altmiu,damagemodfactor,incrementalemissions,liability*(t(2)-t(1))*((plusemissionsplusnonCO2.Consumption-baseline.Consumption)-(plusemissionsplusnonCO2.ROI.Consumption-baseline.ROI.Consumption)));
	end	

	
	
	deltaemissionsplusnonCO2 = (plusemissionsplusnonCO2.Welfare-baseline.Welfare)/sum(incrementalemissions);

	SCCout.SCGHG_consumptiondenominated = (12/44)*-deltaemissionsplusnonCO2./deltaconsumption*1e3;


	Edeltaemissions = (plusemissions.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalemissions);

	EdeltaemissionsplusnonCO2 = (plusemissionsplusnonCO2.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalemissions);

	Edeltaconsumption = (plusconsumption.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalconsumption);

	Edeltacapital = (pluscapital.ExpectedWelfare-baseline.ExpectedWelfare)/sum(incrementalcapital);

	Econvemissions = (plusemissions.EquivalentMaterialExpectedWelfare-baseline.EquivalentMaterialExpectedWelfare)/sum(incrementalemissions);
	Econvconsumption = (plusconsumption.EquivalentMaterialExpectedWelfare-baseline.EquivalentMaterialExpectedWelfare)/sum(incrementalemissions);
	SCCout.ESCC_consumptiondenominated = (12/44)*-Edeltaemissions./Edeltaconsumption*1e3;

	SCCout.ESCC_investmentdenominated = (12/44)*-Edeltaemissions./Edeltacapital*1e3;

	SCCout.ESCC_materialconsumptionequiv = (12/44)*-Econvemissions./Econvconsumption*1e3;;
	SCCout.ESCGHG_consumptiondenominated = (12/44)*-EdeltaemissionsplusnonCO2./Edeltaconsumption*1e3;


	Edeltaabatement = (plusabatement.ExpectedWelfare-baseline.ExpectedWelfare)./mean(sum(incrementalabatement,2));
	SCCout.EMAC_consumptiondenominated = (12/44)*-Edeltaabatement./Edeltaconsumption*1e3;
	SCCout.EMAC_investmentdenominated = (12/44)*-Edeltaabatement./Edeltacapital*1e3;
	
	SCCout.deltaemissions=deltaemissions;
	SCCout.deltaconsumption=deltaconsumption;
	SCCout.deltacapital=deltacapital;
	SCCout.deltaabatement=deltaabatement;
	SCCout.Edeltaemissions=Edeltaemissions;
	SCCout.Edeltaconsumption=Edeltaconsumption;
	SCCout.Edeltacapital=Edeltacapital;
	SCCout.Edeltaabatement=Edeltaabatement;
	
	
	if nargout > 1
		scenarios.baseline = baseline;
		scenarios.plusemissions = plusemissions;
		scenarios.plusconsumption = plusconsumption;
		scenarios.pluscapital = pluscapital;
		scenarios.plusabatement = plusabatement;
	end

end
