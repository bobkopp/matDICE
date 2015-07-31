function [ExpectedWelfare,scenario] = DICEEconomicModel(p,SavingsRate,miu,damagemodfactor,incrementalemissions,incrementalconsumption,incrementalcapital,incrementalforcing,incrementalmiu)

%[ExpectedWelfare,scenario] = DICEEconomicModel(p,SavingsRate,miu,[damagemodfactor],[incrementalemissions],[incrementalconsumption],[incrementalcapital],[incrementalforcing],[incrementalmiu])
%
%
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Mar 4 17:54:50 EST 2013

	Nscenarios = length(p.T2xCO2);

	t=p.t;
	lam = p.lam(:);
	Population = repmat(p.L,Nscenarios,1);
	Nt = length(t);
	
	SavingsRate=SavingsRate(:)';
	if size(miu,1)==1
		miu=repmat(miu(:)',Nscenarios,1);
	end

	defval('damagemodfactor',1);
	defval('incrementalemissions',zeros(size(t)));
	defval('incrementalconsumption',zeros(size(t)));
	defval('incrementalcapital',zeros(size(t)));
	defval('incrementalforcing',zeros(size(t)));
	defval('incrementalmiu',zeros(size(t)));

	if p.doROI
		ROI_Population = bsxfun(@times,Population,p.ROI_fpop);
	else
		ROI_Population = Population;
	end

	elasmu = p.elasmu;
	mincons=p.mincons;
	if elasmu == 1
		calcutility = @(ConsumptionPC) real ( log(max(mincons/1000,ConsumptionPC)) );
		elasmu = 1.0001;
	else
		calcutility = @(ConsumptionPC) real ( ((max(mincons/1000,ConsumptionPC)).^(1-elasmu))./(1-elasmu) );
	end

	b = p.b^((t(2)-t(1))/10);
	c1 = p.c1;
	c3 = p.c3;
	c4 = p.c4;
	gama = p.gama;
	expcost2 = p.expcost2;
	FCO22x = p.FCO22x;
	matPI = p.matPI;
	plantlife=p.plantlife;
	earlyretcost = p.earlyretcost;
	rr=p.rr;
	dk=p.dk;
	cost1 = p.cost1;
	popfeedback=p.popfeedback;
	allowrecoveryfromcrash = p.allowrecoveryfromcrash;
	limmiu = p.limmiu;
	Tstoch = p.Tstoch;
	if Tstoch>0
		Tstoch_seeds = p.Tstoch_seeds;
	end
	
	
	doROI = p.doROI;
	if doROI
		ROI_fcons = p.ROI_fcons;
		ROI_fmiu = p.ROI_fmiu;
		ROI_femit = p.ROI_femit;
		ROI_fliability = p.ROI_fliability;
		ROI_damfactor = p.ROI_damfactor;
		liabilitydamagesref = p.liabilitydamagesref;
		ROI_limmiu = p.ROI_limmiu;
		ROI_anchor = p.ROI_anchor;
		ROI_Tstoch = p.ROI_Tstoch;
		if ROI_Tstoch>0
			ROI_Tstoch_seeds = p.ROI_Tstoch_seeds;
		end
	end

	ecoshare = p.ecoshare;
	
	partfract = p.partfract;
	minoutpc = p.minoutpc;
	
	fosslim = p.fosslim;
	al = p.al;
	etree = p.etree;
	
	consttcre = p.consttcre;
	
	sigma=p.sigma;
	dam_futility = p.dam_futility;
	dam_fcapital = p.dam_fcapital;
	abate_futility = p.abate_futility;
	
	forcoth=p.forcoth;
	forcothmiufactor=p.forcothmiufactor;
	forcothlimmiu=p.forcothlimmiu;
	liabilitycap = p.liabilitycap;
	liabilityfloor = p.liabilityfloor;
	
	miu = bsxfun(@plus,miu,incrementalmiu);
	
	
	calcnextmassdist = @(mass0,emissions) mass0 * b + [emissions zeros(Nscenarios,2)]; 
	if p.use_new_climate
		dt = t(2)-t(1);
		sbconst = p.sbconst;
		Cp_mixed = p.Cp_mixed;
		z_mixed = p.z_mixed;
		z_deep = p.z_deep;
		ocean_diffusivity = p.ocean_diffusivity;
		FCO22x = p.FCO22x;
		Ndeepboxes = p.Ndeepboxes;
		calcnexttemperature = @(Atm,Ocean,Forcing)  DCM(Forcing,Atm,Ocean,FCO22x./lam,ocean_diffusivity,sbconst,Cp_mixed,z_mixed,z_deep,Ndeepboxes,dt);
	else
		Ndeepboxes = 1;
		if consttcre == 0
			calcnexttemperature = @(Atm,Ocean,Forcing) [Atm+c1*(Forcing-lam.*Atm-c3*(Atm-Ocean)) Ocean+c4*(Atm-Ocean)];
		else
			calcnexttemperature = @(Atm,Ocean,Emissions) [Atm+consttcre*Emissions Ocean+c4*(Atm-Ocean)];
		end
	end
	calcoutput_gross = @(TFP,Kapital,Labor) TFP .* Labor.^(1-gama) .* Kapital.^(gama);
	calcabatecost = @(PartFraction,CostFactor,Mitigation) (PartFraction^(1-expcost2))*(bsxfun(@times,CostFactor,Mitigation.^expcost2));
	calcco2forcing = @(AtmCO2) FCO22x*log2((AtmCO2+1e-9)/matPI);
	calcearlyretcost = @(DeltaMitigation,sigma,dt) (abs(DeltaMitigation)>dt./plantlife).*0.5*(max(DeltaMitigation-dt./plantlife).^2).*plantlife.*sigma.*earlyretcost;
	if isprop(p,'calcdamagesstr')
		eval(['calcdamages = ' p.calcdamagesstr ';']);
	else
		calcdamages = p.calcdamages;
	end

	Tatm(:,1) = ones(Nscenarios,1)*p.Tatm0;
	preTatm = Tatm(:,1) - (t(2)-t(1))*p.Trate0;
	if length(p.Tocean0)==1
		Tocean = repmat(reshape(p.Tocean0,1,1,[]),[Nscenarios,1,Ndeepboxes]);
	else
		Tocean = repmat(reshape(p.Tocean0,1,1,[]),[Nscenarios,1,1]);
	end
	MassDist(1,1:3,1:Nscenarios) = repmat([p.mat2005 p.mu2005 p.ml2005],[1 1 Nscenarios]);
	ForcingNonCO2(:,1) = repmat(forcoth(1),Nscenarios,1)+incrementalforcing(:,1);
	Forcing(:,1) = squeeze(calcco2forcing(MassDist(1,1,:)))+ForcingNonCO2(:,1);
	CumulativeEmissions(:,1) = zeros(Nscenarios,1);
	dTatm(:,1) = ones(Nscenarios,1)*p.Trate0;

	Capital(:,1) = ones(Nscenarios,1)*p.k0;
	CapitalClimateDamages(:,1) = zeros(size(Capital(:,1)));
	Capital(:,1) = Capital(:,1) + incrementalcapital(1);
	EarlyRetirementCost(:,1) = zeros(Nscenarios,1);
	NotCrashed(:,1) = ones(Nscenarios,1);
	
for i=1:length(t)
	
		if i>1
			if Nscenarios > 1
				MassDist(i,1:3,1:Nscenarios) = reshape(calcnextmassdist(squeeze(MassDist(i-1,:,:))',Emissions(:,i-1))',1,3,Nscenarios);
				MassDist(i,:,:) = MassDist(i,:,:).*(MassDist(i,:,:)>0);
			else
				MassDist(i,1:3,1) = calcnextmassdist(MassDist(i-1,:),Emissions(i-1));
				MassDist(i,1:3,1) = MassDist(i,:).*(MassDist(i,:)>0);
			end
			forcothmiu(:,i) = miu(:,i-1)*forcothmiufactor;
			forcothmiu(:,i) = (forcothmiu(:,i)<=forcothlimmiu).*forcothmiu(:,i) + forcothlimmiu.*(forcothmiu(:,i)>forcothlimmiu);
			ForcingNonCO2(:,i)=real((forcoth(i)*(1-forcothmiu(:,i)))+incrementalforcing(:,i));
			Forcing(:,i) = real(calcco2forcing(squeeze(MassDist(i,1,:)))+ForcingNonCO2(:,i));
			if consttcre == 0
				temps = real(calcnexttemperature(Tatm(:,i-1),Tocean(:,i-1,:),.5*(Forcing(:,i)+Forcing(:,i-1))));
				CumulativeEmissions(:,i) = CumulativeEmissions(:,i-1) + Emissions(:,i-1);
			else
				temps = real(calcnexttemperature(Tatm(:,i-1),Tocean(:,i-1,:),Emissions(:,i-1)));
				CumulativeEmissions(:,i) = CumulativeEmissions(:,i-1) + Emissions(:,i-1);
			end

			% hard upper limit to carbon
			miu(:,i) = miu(:,i) .* (CumulativeEmissions(:,i)<fosslim) + (CumulativeEmissions(:,i)>=fosslim);
			EarlyRetirementCost(:,i) = calcearlyretcost((miu(:,i)-miu(:,i-1))/partfract(i),sigma(i),t(2)-t(1));

			Tatm(:,i) = temps(:,1); Tocean(:,i,:) = permute(temps(:,2:end),[1 3 2]);
			if Tstoch>0
				Tatm(:,i) = Tatm(:,i) + Tstoch_seeds(:,i);
			end
			dTatm(:,i) = (Tatm(:,i)-Tatm(:,i-1))/(t(i)-t(i-1));
			
			Capital(:,i) = Capital(:,i-1)*(1-dk).^(t(i)-t(i-1)) + Investment(:,i-1)*(t(i)-t(i-1));
			CapitalClimateDamages(:,i) = 1 - (1 - damagemodfactor*calcdamages(abs(Tatm(:,i)),t(i),dTatm(:,i),Output_Gross(:,i-1),[preTatm Tatm],Population(i))).^(dam_fcapital./gama);
			Capital(:,i) = Capital(:,i) .* (1 - CapitalClimateDamages(:,i));
			Capital(:,i) = Capital(:,i) + incrementalcapital(i);
		end			
			
			
	
		Output_Gross(:,i) = calcoutput_gross(al(:,i),Capital(:,i),Population(:,i));
		
		FossilEmissions_Gross(:,i) = (t(2)-t(1))*sigma(i).*Output_Gross(:,i);
		FossilEmissions(:,i) = FossilEmissions_Gross(:,i).*(1-miu(:,i));
		LandEmissions_Gross(:,i) = repmat(etree(i)*(t(2)-t(1))/10,Nscenarios,1);
		LandEmissions(:,i) = LandEmissions_Gross(:,i) .*(1-miu(:,i));
		Emissions(:,i) = FossilEmissions(:,i) + LandEmissions(:,i) + incrementalemissions(i);  

		ClimateDamages(:,i) = 1 - (1 - damagemodfactor*calcdamages(abs(Tatm(:,i)),t(i),dTatm(:,i),Output_Gross(:,i),[preTatm Tatm],Population(i))).^(1-dam_fcapital-dam_futility);
		AbatementCost(:,i) = calcabatecost(partfract(i),cost1(i),miu(:,i)) + EarlyRetirementCost(:,i);
		
		AbatementCostOutput(:,i) = 1 - (1 - AbatementCost(:,i)).^(1-abate_futility);
		AbatementCostUtility(:,i) = 1 - (1 - AbatementCost(:,i)).^(abate_futility);

		Output(:,i) = Output_Gross(:,i) .* (1-ClimateDamages(:,i)-AbatementCostOutput(:,i));
		minOutput = minoutpc*1e-6*Population(:,i);
		if popfeedback
			Population(:,i) = Population(:,i).*(Output(:,i)>minOutput) + (Output(:,i)/(minoutpc*1e-6)).*(Output(:,i)<=minOutput);
		else
			if (i>1) && (~allowrecoveryfromcrash)
				Output(:,i)=Output(:,i).*NotCrashed(:,i-1); 
			end
			Output(:,i) = Output(:,i).*(Output(:,i)>minOutput) + minOutput.*(Output(:,i)<=minOutput); % ensure output doesn't fall below $730/person
		end
		NotCrashed(:,i) = Output(:,i)>minOutput;
		Investment(:,i) = Output(:,i) * SavingsRate(i);
		
		EcoClimateDamages(:,i) = 1 - (1 - damagemodfactor*calcdamages(abs(Tatm(:,i)),t(i),dTatm(:,i),Output_Gross(:,i),[preTatm Tatm],Population(i))).^(1-dam_futility);
		EffectiveConsumptionDamages(:,i) = 1 - (1 - damagemodfactor*calcdamages(abs(Tatm(:,i)),t(i),dTatm(:,i),Output_Gross(:,i),[preTatm Tatm],Population(i))).^(dam_futility);
	end
	
	Consumption = Output.*repmat((1-SavingsRate),Nscenarios,1);
	if size(incrementalconsumption,1)==1
		Consumption = Consumption + repmat(incrementalconsumption/(t(2)-t(1)),Nscenarios,1);
	else
		Consumption = Consumption + incrementalconsumption;
	end
	ConsumptionPerCapita = 1000*Consumption./Population;
	EcoConsumptionPerCapita = repmat(ConsumptionPerCapita(:,1),1,size(EcoClimateDamages,2)).*(1-EcoClimateDamages);

	if doROI
		ROI_minOutput = minoutpc*1e-6*ROI_Population;
		ROI_Abatement = max(0,bsxfun(@plus,min(ROI_limmiu,bsxfun(@plus,ROI_anchor(2,:),bsxfun(@times,bsxfun(@minus,miu,bsxfun(@plus,ROI_anchor(1,:),incrementalmiu)),ROI_fmiu./ROI_femit))),bsxfun(@rdivide,incrementalmiu,ROI_femit)));
%		ROI_AbatementCost = bsxfun(@times,AbatementCost,(bsxfun(@times,ROI_Abatement,ROI_femit./ROI_fcons)./(miu+eps)));
%		ROI_AbatementCost = bsxfun(@times,AbatementCost,ROI_fmiu./ROI_fcons);
		ROI_AbatementCost = calcabatecost(1,cost1,ROI_Abatement) ;
		ROI_AbatementCostOutput = 1 - (1 - ROI_AbatementCost).^(1-abate_futility);
		ROI_AbatementCostUtility = 1 - (1 - ROI_AbatementCost).^(abate_futility);

		ROI_Emissions_Gross = bsxfun(@times,FossilEmissions_Gross+LandEmissions_Gross,ROI_femit);
		ROI_Emissions = bsxfun(@plus,bsxfun(@times,ROI_Emissions_Gross,1-ROI_Abatement),incrementalemissions);
		%ROI_OutputGross = bsxfun(@times,Output+Output_Gross.*AbatementCost,ROI_fcons);
		ROI_OutputGross = bsxfun(@times,Output_Gross,ROI_fcons);
		if ROI_Tstoch == 0
			ROI_Tatm = Tatm;
			ROI_dTatm = dTatm;
			ROI_ClimateDamages = ClimateDamages * ROI_damfactor;
		else
			ROI_Tatm = Tatm + ROI_Tstoch_seeds;
			ROI_dTatm = [dTatm(:,1) diff(ROI_Tatm,[],2)];
			for i=1:length(t)
					ROI_ClimateDamages(:,i) = 1 - (1 - damagemodfactor*calcdamages(abs(ROI_Tatm(:,i)),t(i),ROI_dTatm(:,i),Output_Gross(:,i),[preTatm ROI_Tatm(:,1:i)],Population(i))).^(1-dam_fcapital-dam_futility);
			end
			ROI_ClimateDamages = ROI_damfactor * ROI_ClimateDamages;
		end
		ROI_Output = max(ROI_OutputGross.*(1-ROI_AbatementCostOutput-ROI_ClimateDamages),ROI_minOutput);

		ROI_Liability = max(liabilityfloor,bsxfun(@times,ROI_fliability,bsxfun(@times,Output_Gross,1-ROI_fcons).*(1-(1-ClimateDamages).*(1-CapitalClimateDamages).*(1-EffectiveConsumptionDamages))));
		%ROI_Liability = -bsxfun(@times,ROI_fliability,bsxfun(@times,Output,1-ROI_fcons));

		%ROI_Liability = min(ROI_Liability,liabilitycap*ROI_OutputGross);
		%ROI_Liability = max(ROI_Liability,-liabilitycap*ROI_OutputGross);
		ROI_Consumption = (ROI_Output) .* repmat((1-SavingsRate),Nscenarios,1);
		if size(incrementalconsumption,1)==1
			ROI_Consumption = ROI_Consumption + repmat(incrementalconsumption/(t(2)-t(1)),Nscenarios,1);
		else
			ROI_Consumption = ROI_Consumption + incrementalconsumption;
		end
		ROI_ConsumptionPerCapita = 1000*(ROI_Consumption)./ROI_Population;

		ROI_EcoConsumptionPerCapita = repmat(ROI_ConsumptionPerCapita(:,1),1,size(EcoClimateDamages,2)).*(1-EcoClimateDamages);
		ROI_LiabilityUtilityPerCapita = -real( max(mincons/1000,ROI_ConsumptionPerCapita).^(-elasmu) .* ( 1000 * ROI_Liability ./ ROI_Population .* repmat((1-SavingsRate),Nscenarios,1)) );
	else
		ROI_ConsumptionPerCapita = ConsumptionPerCapita;
		ROI_EcoConsumptionPerCapita = EcoConsumptionPerCapita;
		ROI_LiabilityUtilityPerCapita = zeros(size(ConsumptionPerCapita));
		ROI_AbatementCostUtility = AbatementCostUtility;
	end

	ppmCO2=squeeze(MassDist(:,1,:))' /2.13;

	if ecoshare==0

		EffectiveConsumptionPerCapita = max(ConsumptionPerCapita .* (1 - EffectiveConsumptionDamages - AbatementCostUtility),mincons/1e3);
		EquivalentMaterialConsumptionPerCapita = EffectiveConsumptionPerCapita;
		ROI_EffectiveConsumptionPerCapita = max(ROI_ConsumptionPerCapita .* (1 - EffectiveConsumptionDamages - ROI_AbatementCostUtility),mincons/1e3);
		ROI_EquivalentMaterialConsumptionPerCapita = ROI_EffectiveConsumptionPerCapita;

		UtilityPC = calcutility(ROI_EffectiveConsumptionPerCapita);
		EquivalentMaterialUtilityPC = UtilityPC;

		UtilityDiscountedTotal = real((UtilityPC + ROI_LiabilityUtilityPerCapita).* ROI_Population .* repmat(rr,Nscenarios,1));
		Welfare = sum(UtilityDiscountedTotal,2);
		ExpectedWelfare = mean(Welfare);

		EquivalentMaterialWelfare = Welfare;
		EquivalentMaterialExpectedWelfare = ExpectedWelfare;

	else

		EffectiveConsumptionPerCapita = ((1-ecosh).*ConsumptionPerCapita.^(1-1./ecoel) + ecosh .* EcoConsumptionPerCapita.^(1-1./ecoel)).^(ecoel./(ecoel-1));
		EffectiveConsumptionPerCapita = max(mincons/1e3,EffectiveConsumptionPerCapita .* (1 - EffectiveConsumptionDamages - AbatementCostUtility));
		EquivalentMaterialConsumptionPerCapita = ( ( ( (1-ecosh).*ConsumptionPerCapita.^(1-1./ecoel) + ecosh .* EcoConsumptionPerCapita.^(1-1./ecoel) - ecosh .* repmat(EcoConsumptionPerCapita(:,1).^(1-1./ecoel(:,1)),1,length(t)) ) ./ (1-ecosh)) .^ (ecoel./(ecoel-1)) ) .* (1 - EffectiveConsumptionDamages);
		ROI_EffectiveConsumptionPerCapita = ((1-ecosh).*ROI_ConsumptionPerCapita.^(1-1./ecoel) + ecosh .* ROI_EcoConsumptionPerCapita.^(1-1./ecoel)).^(ecoel./(ecoel-1));
		ROI_EffectiveConsumptionPerCapita = max(mincons/1e3,EffectiveConsumptionPerCapita .* (1 - ROI_EffectiveConsumptionDamages - ROI_AbatementCostUtility));
		ROI_EquivalentMaterialConsumptionPerCapita = ( ( ( (1-ecosh).*ROI_ConsumptionPerCapita.^(1-1./ecoel) + ecosh .* ROI_EcoConsumptionPerCapita.^(1-1./ecoel) - ecosh .* repmat(ROI_EcoConsumptionPerCapita(:,1).^(1-1./ecoel(:,1)),1,length(t)) ) ./ (1-ecosh)) .^ (ecoel./(ecoel-1)) ) .* (1 - ROI_EffectiveConsumptionDamages);

		UtilityPC = calcutility(ROI_EffectiveConsumptionPerCapita);
		EquivalentMaterialUtilityPC = calcutility(ROI_EquivalentMaterialConsumptionPerCapita);

		UtilityDiscountedTotal = real((UtilityPC + ROI_LiabilityUtilityPerCapita).* ROI_Population .* repmat(rr,Nscenarios,1));
		Welfare = sum(UtilityDiscountedTotal,2);
		ExpectedWelfare = mean(Welfare);

		EquivalentMaterialWelfare = sum(real(EquivalentMaterialUtilityPC .* ROI_Population .* repmat(rr,Nscenarios,1)),2);
		EquivalentMaterialExpectedWelfare = mean(EquivalentMaterialWelfare);

	end

	% calculate discount rates
	Growth = ROI_ConsumptionPerCapita./repmat(ROI_ConsumptionPerCapita(:,1),1,size(ConsumptionPerCapita,2));
	EffectiveGrowth = ROI_EffectiveConsumptionPerCapita./repmat(ROI_EffectiveConsumptionPerCapita(:,1),1,size(ConsumptionPerCapita,2));
	DiscountFactor = repmat(rr,Nscenarios,1).*(Growth.^-elasmu);
	EffectiveDiscountFactor = repmat(rr,Nscenarios,1).*(EffectiveGrowth.^-elasmu);
	DiscountRate = bsxfun(@power,DiscountFactor(:,1:end-1)./DiscountFactor(:,2:end),1./diff(t))-1;
	EffectiveDiscountRate = bsxfun(@power,EffectiveDiscountFactor(:,1:end-1)./EffectiveDiscountFactor(:,2:end),1./diff(t))-1;

	Emissions_Gross = bsxfun(@plus,FossilEmissions_Gross + LandEmissions_Gross,incrementalemissions);

	if nargout>1

		scenario.Population = Population;
		scenario.ExpectedWelfare = ExpectedWelfare;
		scenario.Welfare = Welfare;
		scenario.abatement = miu;
		scenario.abatementnonco2 = forcothmiu;
		scenario.Consumption = Consumption;
		scenario.ConsumptionPerCapita = ConsumptionPerCapita;
		scenario.EcoConsumptionPerCapita = EcoConsumptionPerCapita;
		scenario.EffectiveConsumptionPerCapita = EffectiveConsumptionPerCapita;
		scenario.UtilityPerCapita = UtilityPC;
		scenario.UtilityDiscountedTotal = UtilityDiscountedTotal;
		scenario.Output_Gross = Output_Gross;
		scenario.Output = Output;
		scenario.Investment = Investment;
		scenario.Capital = Capital;
		scenario.ClimateDamages = ClimateDamages;
		scenario.EcoClimateDamages = EcoClimateDamages;
		scenario.CapitalClimateDamages = CapitalClimateDamages;
		scenario.EffectiveConsumptionDamages = EffectiveConsumptionDamages;
		scenario.AbatementCost = AbatementCost;
		scenario.EarlyRetirementCost = EarlyRetirementCost;
		scenario.FossilEmissions = FossilEmissions;
		scenario.Emissions = Emissions;
		scenario.Emissions_Gross = Emissions_Gross;
		scenario.Tatm = Tatm;
		scenario.Tocean = Tocean;
		scenario.ppmCO2 = ppmCO2;
		scenario.Forcing = Forcing;
		scenario.ForcingNonCO2 = ForcingNonCO2;
	%	scenario.UtilityScalePoints = [utilscalptsy(:) utilscalptsx(:)];
		scenario.EquivalentMaterialConsumptionPerCapita = EquivalentMaterialConsumptionPerCapita;
		scenario.EquivalentMaterialUtilityPerCapita = EquivalentMaterialUtilityPC;
		scenario.EquivalentMaterialWelfare = EquivalentMaterialWelfare;
		scenario.EquivalentMaterialExpectedWelfare = EquivalentMaterialExpectedWelfare;
		scenario.DiscountFactor = DiscountFactor;
		scenario.EffectiveDiscountFactor = EffectiveDiscountFactor;
		scenario.DiscountRate = DiscountRate;
		scenario.EffectiveDiscountRate = EffectiveDiscountRate;
		if doROI
			scenario.ROI.Population = ROI_Population;
			scenario.ROI.Consumption = ROI_Consumption;
			scenario.ROI.Abatement = ROI_Abatement;
			scenario.ROI.AbatementCost = ROI_AbatementCost;
			scenario.ROI.Emissions = ROI_Emissions;
			scenario.ROI.Emissions_Gross = bsxfun(@times,Emissions_Gross,ROI_femit);
			scenario.ROI.Liability = ROI_Liability;
			scenario.ROI.Output = ROI_Output;
			scenario.ROI.Output_Gross = ROI_OutputGross;
			scenario.ROI.ClimateDamages = ROI_ClimateDamages;
			scenario.ROI.EffectiveConsumptionPerCapita = ROI_ConsumptionPerCapita;
			scenario.ROI.LiabilityUtilityPerCapita = ROI_LiabilityUtilityPerCapita;
			scenario.ROI.Tatm = ROI_Tatm;
		else
			scenario.ROI = 0;
		end

		%scenario.CarbonEmittingEnergyStock = bsxfun(@times,sigma,Output_Gross.*max(0,1-miu));
		%scenario.CarbonFreeEnergyStock = bsxfun(@times,sigma,Output_Gross.*min(1,miu));
		%scenario.CarbonNegativeEnergyStock = bsxfun(@times,sigma,Output_Gross.*max(0,miu-1));
		%scenario.NewBuildsCarbonEmitting(:,2:Nt) = scenario.CarbonEmittingEnergyStock(:,2:Nt) - %scenario.CarbonEmittingEnergyStock(:,1:Nt-1).*(t(2)-t(1))./plantlife;
		%scenario.NewBuildsCarbonFree(:,2:Nt) = scenario.CarbonFreeEnergyStock(:,2:Nt) - %scenario.CarbonFreeEnergyStock(:,1:Nt-1).*(t(2)-t(1))./plantlife;
		%scenario.NewBuildsCarbonNegative(:,2:Nt) = scenario.CarbonNegativeEnergyStock(:,2:Nt) - %scenario.CarbonNegativeEnergyStock(:,1:Nt-1).*(t(2)-t(1))./plantlife;
		
	end
	
end

function newTemps = DCM(Forcing,Atm,Ocean,cses,ocean_diffusivity,sbconst,Cp_mixed,z_mixed,z_deep,Ndeepboxes,dt)
		y = DiffusiveClimateModel(Forcing,[Atm permute(Ocean,[1 3 2])],cses,ocean_diffusivity,sbconst,Cp_mixed,z_mixed,z_deep,Ndeepboxes,dt);
	newTemps=permute(y(end,:,:),[3 2 1]);
end
