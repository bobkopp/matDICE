
function [Welfare,Consumption,ConsumptionPerCapita,EcoConsumptionPerCapita,InstantaneousUtilityPC,Output_Gross,Output,Investment,Capital,ClimateDamages,AbatementCost,Emissions,CumulativeEmissions,Tatm,Tocean,ppmCO2,Forcing] = RICEEconomicModel(p,SavingsRate,miu,varargin)

	t=p.t;
	nreg = length(p.q0);

	calcnextmassdist = @(mass0,emissions) mass0 * p.b + [emissions 0 0];
	calcnexttemperature = @(Atm,Ocean,Forcing) [Atm+p.c1*(Forcing-p.lam*Atm-p.c3*(Atm-Ocean)) Ocean+p.c4*(Atm-Ocean)];
	calcco2forcing = @(AtmCO2) p.FCO22x*log2((AtmCO2+1e-9)/p.matPI);
	
	calcoutput_gross = @(TFP,Kapital,Labor) TFP .* Labor.^(1-p.gama) .* Kapital.^(p.gama);
	calcdamfunc = @(Temperature) p.aa1*Temperature+p.aa2*Temperature.^p.aa3;
	calcdamages = @(Temperature) calcdamfunc(Temperature) ./ (1+calcdamfunc(Temperature));
	calcabatecost = @(CostFactor,Mitigation) (CostFactor.*Mitigation.^p.expcost2);

	damagemodfactor = 1;
	incrementalemissions = zeros(nreg,length(t));
	incrementalconsumption = zeros(nreg,length(t));
	incrementalcapital = zeros(nreg,length(t));
	if length(varargin)>0
		damagemodfactor = varargin{1};
	end
	if length(varargin)>1
		incrementalemissions = varargin{2};
	end
	if length(varargin)>2
		incrementalconsumption = varargin{3};
	end
	if length(varargin)>3
		incrementalcapital = varargin{4};
	end
	
	Tatm(1) = p.Tatm0;
	Tocean(1) = p.Tocean0;
	MassDist(1,:) = [p.mat2005 p.mu2005 p.ml2005];
	Forcing(1) = calcco2forcing(MassDist(1,1))+p.forcoth(1);
	CumulativeEmissions(1) = 0;

	Capital(:,1) = p.k0' + incrementalcapital(:,1);
	
	for i=1:length(t)
	
		if i>1
			MassDist(i,:) = calcnextmassdist(MassDist(i-1,:),Emissions(i-1));
			Forcing(i) = calcco2forcing(MassDist(i,1))+p.forcoth(i);
			temps = calcnexttemperature(Tatm(i-1),Tocean(i-1),.5*(Forcing(i)+Forcing(i-1)));
			CumulativeEmissions(i) = CumulativeEmissions(i-1) + sum(Emissions(:,i-1));

			Tatm(i) = temps(1); Tocean(i) = temps(2);
			
			Capital(:,i) = Capital(:,i-1)*(1-p.dk).^(t(i)-t(i-1)) + Investment(:,i-1)*(t(i)-t(i-1)) + incrementalcapital(:,i);
		end				

		Output_Gross(:,i) = calcoutput_gross(p.al(:,i),Capital(:,i),p.L(:,i));
		% hard upper limit to carbon
		miu(:,i) = miu(:,i) * (CumulativeEmissions(i)<p.fosslim) + ones(size(miu(:,i))) * (CumulativeEmissions(i)>=p.fosslim);
		
		Emissions(:,i) = (t(2)-t(1))*p.sigma(:,i).*(1-miu(:,i)).*Output_Gross(:,i) + p.etree(:,i) + incrementalemissions(:,i);
		ClimateDamages(:,i) = damagemodfactor.*calcdamages(Tatm(i));
		AbatementCost(:,i) = calcabatecost(p.cost1(:,i),miu(:,i));
		Output(:,i) = Output_Gross(:,i) .* (1-ClimateDamages(:,i)-AbatementCost(:,i));
		Investment(:,i) = Output(:,i) .* SavingsRate(:,i);
		
	end
	
	Consumption = Output.*(1-SavingsRate) + incrementalconsumption; 
	ConsumptionPerCapita = 1000*Consumption./p.L;
	
	EcoAbundance = 1./(1+p.ecodamagec1*Tatm.^p.ecodamagec2);
	EcoConsumptionPerCapita = ConsumptionPerCapita(1)*EcoAbundance/EcoAbundance(1);
	
	ppmCO2=MassDist(:,1)' /2.13;

	InstantaneousUtilityPC0 = 1 + (ConsumptionPerCapita.^(1-p.elasmu))./(1-p.elasmu);
	InstantaneousUtilityPC = InstantaneousUtilityPC0;
	%InstantaneousUtilityPC = (( (1-p.ecoshare)*ConsumptionPerCapita.^(1-1./p.ecoelasticity) + p.ecoshare * EcoConsumptionPerCapita.^(1-1./p.ecoelasticity)).^((1-p.elasmu)*p.ecoelasticity/(p.ecoelasticity-1)))./(1-p.elasmu);
	%InstantaneousUtilityPC = InstantaneousUtilityPC * InstantaneousUtilityPC0(1)/InstantaneousUtilityPC(1);
	
	InstantaneousUtilityTotalPeriod = bsxfun(@times,InstantaneousUtilityPC .* p.L, p.rr);
	InstantaneousUtilityWorld = sum(InstantaneousUtilityTotalPeriod,1);
	Trajectory = .5*(sum(diff(t).*(InstantaneousUtilityWorld(1:end-1)+InstantaneousUtilityWorld(2:end))));
	%Trajectory = 10*sum(InstantaneousUtilityTotalPeriod);
	Welfare = real(Trajectory*p.scale1 + p.scale2);

end