if matlabpool('size')==0
	matlabpool('open',feature('numCores'));
end


% Load MiniCAM based scenario

	defp=DICEParameters;
	years = 2000:10:2300;
	GDP = [36.1 47.4 60.8 78.9 100.6 125.7 160.6 201.8 249.3 309.4 369.5 437.8 514.6 600.2 694.5 797.3 908.2 1026.3 1150.6 1279.8 1412.4 1543.4 1670.0 1789.3 1898.3 1994.2 2074.4 2136.7 2179.2 2200.8 2200.8];
	Population = [6.04 6.80 7.53 8.07 8.50 8.81 8.98 9.04 9.01 8.84 8.66 8.51 8.38 8.27 8.17 8.09 8.03 7.98 7.95 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93 7.93]*1000;
	IndustrialCO2 = [26.53 31.81 38.01 45.14 51.74 57.80 62.75 67.78 72.89 76.67 80.46 83.77 86.53 88.68 90.17 90.96 91.03 90.39 89.05 87.04 84.40 81.04 77.06 72.55 67.63 62.43 57.06 51.65 46.28 41.07 36.09] * 12/44;
	LandCO2 = [3.82 2.80 1.93 1.84 1.77 1.73 1.62 1.32 0.85 0.73 0.62 0.56 0.50 0.43 0.37 0.31 0.25 0.19 0.12 0.06 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
	nonCO2forcings = [0.67 0.75 0.85 0.94 1.03 1.13 1.19 1.25 1.31 1.34 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41 1.41];
	
	[Y,L,E,K,eland,forcoth,al,sigma]=DICEFitRefScenario(years,GDP,Population,IndustrialCO2,LandCO2,nonCO2forcings);
	
	specMCam={'q0',Y(1),'k0',K(1),'a0',al(1),'sig0',sigma(1),'pop0',L(1),'L',L,'al',al,'sigma',sigma,'etree',eland,'forcoth',forcoth};

% discount parameters

	Ref(1) = matDICE(specMCam{:},'calcdamages','0*Temperature');
	
	ConsPCgrowth = (Ref(1).ConsumptionPerCapita(12)/Ref(1).ConsumptionPerCapita(2))^(1/(Ref(1).times(12)-Ref(1).times(2)));
	targetdiscount = 0.03;
	etas=[0 1 1.4 2];
	rhos=(1+targetdiscount)./(ConsPCgrowth.^etas)-1;

	for i=1:length(etas)
		spec_disc{i} = {'prstp',rhos(i),'elasmu',etas(i)};
	end
	
	minoutpcs = [500 125 2000];
	for i=1:length(minoutpcs)
		spec_minout{i} = {'minoutpc', minoutpcs(i),'mincons',minoutpcs(i)};
	end
	
% calibration parameters

	Base(1) = matDICE(specMCam{:});
	calib(1) = Base(1).p.calcdamages(2.5);
	calib(2) = .5*calib(1);
	calib(3) = 2*calib(1);
	
	BaseCalibTime = interp1(Base(1).Tatm(1:10),Base(1).times(1:10),2.5);
	BaseCalibRate = 0.1/(interp1(Base(1).Tatm(1:10),Base(1).times(1:10),2.6)-BaseCalibTime);
	BaseCalibGDP = interp1(Base.Tatm(1:15),Base.Output_Gross(1:15),2.5);
	BaseCalibPop = interp1(Base.Tatm(1:15),Base.Population(1:15),2.5);
	
	BaseCalibIncomePC = BaseCalibGDP/BaseCalibPop;

% 2.5 C stabilization reference

	[a,Stabilization] = matDICE(specMCam{:},spec_disc{1}{:},'calcdamages','1./(1+exp(-(Temperature-2.6)*100))','T2xCO2',3,'monotonicAbatement')
	StabAbate = Stabilization.abatement;
	%StabAbate =  [ 0    0.2310    0.2868    0.3579    0.4421    0.5336    0.6304    0.7349    0.8336    0.8754    0.8754    0.8754    0.8754    0.8754    0.8754    0.8881    0.9005    0.9054    0.9086    0.9139    0.9155    0.9155    0.9155    0.9155    0.9155    0.9155    0.9155    0.9155    0.9155    0.9155    0.9155];


% seeds for sampling

	N=1000;
	clear seq;
	for j=1:20
		seq(:,1,j)=1:N;
		for i=2:8
			seq(:,i,j)=randperm(N);
		end
		cor = corrcoef(seq(:,:,j)).^2;
		cors(j) = sum(cor(:));
	end
	[m,mi]=min(cors);
	seeds = seq(:,:,mi)/N - rand(size(seq,1),size(seq,2))/N;
	cses = icdfRoeBaker(seeds(:,1));

	specCS = {'T2xCO2',cses};

clear labels longlabs spec_dam;
labels={}; longlabs = {};
index=0;

% Reference scenarios with uncertain cs

	RefUnc(1) = matDICE(specMCam{:},'calcdamages','0*Temperature','T2xCO2',cses);
	RefStabUnc(1) = matDICE(specMCam{:},'calcdamages','0*Temperature','baseabatement',StabAbate,'T2xCO2',cses);


% 1. DICE Damage Function
index=index+1;
spec_dam{index,1} = {'aa2',Base(1).p.aa2};
spec_dam{index,2} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(2)/calib(1)};
spec_dam{index,3} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(3)/calib(1)};
labels={labels{:},'D'};
longlabs={longlabs{:},'D - DICE'};

% 7. Weitzman exponential
index=index+1;
spec_dam{index,1} = {'calcdamages','1 - exp(-p.aa2 .* Temperature.^p.aa3)'};
spec_dam{index,2} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(2)/calib(1)};
spec_dam{index,3} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(3)/calib(1)};
labels={labels{:},'We'};
longlabs={longlabs{:},'We - Weitzman exponential'};

% 13. Lempert
index=index+1;
spec_dam{index,1} = {'aa4',0.0033,'aa5',3,'aa2',(calib(1) - 0.0033*(10*BaseCalibRate/0.35)^3)/((2.5/3)^2),'aa3',2,'calcdamages','1-1./(1+p.aa2 .* (Temperature/3).^p.aa3 + p.aa4 .* ( (Temperature - mean(TempHist(:,max(size(TempHist,2)-2,1):end),2) ) / 0.35 ) .^ p.aa5)'};
spec_dam{index,2} = {'aa4',0.0033*calib(2)/calib(1),'aa5',3,'aa2',(calib(2) - 0.0033*calib(2)/calib(1)*(10*BaseCalibRate/0.35)^3)/((2.5/3)^2),'aa3',2,'calcdamages','1-1./(1+p.aa2 .* (Temperature/3).^p.aa3 + p.aa4 .* ( (Temperature - mean(TempHist(:,max(size(TempHist,2)-2,1):end),2) ) / 0.35 ) .^ p.aa5)'};
spec_dam{index,3} = {'aa4',0.0033*calib(3)/calib(1),'aa5',3,'aa2',(calib(3) - 0.0033*calib(3)/calib(1)*(10*BaseCalibRate/0.35)^3)/((2.5/3)^2),'aa3',2,'calcdamages','1-1./(1+p.aa2 .* (Temperature/3).^p.aa3 + p.aa4 .* ( (Temperature - mean(TempHist(:,max(size(TempHist,2)-2,1):end),2) ) / 0.35 ) .^ p.aa5)'};
labels={labels{:},'L'};
longlabs={longlabs{:},'L - Lempert et al.'};


% 6. Sterner-Persson Damage Function
index=index+1;
spec_dam{index,1} = {'ecoshare',0.1,'ecoelasticity',0.5};
spec_dam{index,2} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(2)/calib(1)};
spec_dam{index,3} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(3)/calib(1)};
labels={labels{:},'SP'};
longlabs={longlabs{:},'SP - Sterner and Persson'};

% 14. Weitzman (2009)

index=index+1;
aa4 = 1;
spec_dam{index,1} = {'aa2',Base(1).p.aa2,'calibincomepc',BaseCalibIncomePC,'aa4',aa4,'calcdamages','1-1./(1+(Output./Population/p.calibincomepc).^p.aa4.*p.aa2.*abs(Temperature).^p.aa3)','dam_futility',1};
spec_dam{index,2} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(2)/calib(1),};
spec_dam{index,3} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(3)/calib(1)};
labels={labels{:},'Wa'};
longlabs={longlabs{:},'Wa - Weitzman additive'};


% Ad. FUND-like adaptation damage function
index=index+1;
aa4 = -0.4;
spec_dam{index,1} = {'aa2',Base(1).p.aa2,'calibincomepc',BaseCalibIncomePC,'aa4',aa4,'calcdamages','1-1./(1+(Output./Population/p.calibincomepc).^p.aa4.*p.aa2.*abs(Temperature).^p.aa3)'};
spec_dam{index,2} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(2)/calib(1),};
spec_dam{index,3} = {spec_dam{index,1}{:},'aa2',Base(1).p.aa2*calib(3)/calib(1)};
labels={labels{:},'Ad'};
longlabs={longlabs{:},'Ad - Adaptive'};

% 11. Keller damage function
index=index+1;
aa4s = seeds(:,2)*.03;
co2ecrit = @(ClimateSensitivity) exp(-1.52*(log(ClimateSensitivity)-log(3)) + 6.985); 
Tcrit = @(ClimateSensitivity) (log(co2ecrit(ClimateSensitivity)/280)/log(2))*(2.6/3)*ClimateSensitivity;
aa5 = Tcrit(3); % could instead vary this, but I think it's a little weird as expressed in the paper -- really should vary as a function of temperature, not RF
spec_dam{index,1} = {'aa4', aa4s, 'aa5',aa5,'calcdamages', '1-1./(1+p.aa2*Temperature.^p.aa3 + p.aa4 ./ (1 + exp(-(Temperature - p.aa5)*100)))'};
spec_dam{index,2} = {'aa2',defp.aa2*calib(2)/calib(1),spec_dam{index,1}{:},'aa4',aa4s*calib(2)/calib(1)};
spec_dam{index,3} = {'aa2',defp.aa2*calib(3)/calib(1),spec_dam{index,1}{:},'aa4',aa4s*calib(3)/calib(1)};
labels={labels{:},'K'};
longlabs={longlabs{:},'K - Keller et al.'};


% 10. Azar and Lindgren Damage Function (modified)
index=index+1;

aa3s = ones(N,1)*defp.aa3;
aa3s(round(0.97*N):end) = 4;

aa2s = ones(N,1)*defp.aa2;
aa2s(round(0.97*N):end) = defp.aa2 * 2.5.^(defp.aa3)./(2.5+2.5.^4);

aa1s = zeros(N,1);
aa1s(round(0.97*N):end) = aa2s(round(0.97*N):end);

spec_dam{index,1} = {'aa1',aa1s,'aa2',aa2s, 'aa3', aa3s};
spec_dam{index,2} = {'aa1',aa1s,'aa2',aa2s.*(calib(2)/calib(1)), 'aa3', aa3s};
spec_dam{index,3} = {'aa1',aa1s,'aa2',aa2s.*(calib(3)/calib(1)), 'aa3', aa3s};

labels={labels{:},'AL'};
longlabs={longlabs{:},'AL - Azar and Lindgren'};


% 8. Ackerman damage function
index=index+1;
aa3s = icdftri(seeds(:,2),[1 2 5]);
aa2s = repmat(calib,N,1) ./ (2.5.^repmat(aa3s,1,3));
spec_dam{index,1} = {'aa2',aa2s(:,1), 'aa3', aa3s};
spec_dam{index,2} = {'aa2',aa2s(:,2), 'aa3', aa3s};
spec_dam{index,3} = {'aa2',aa2s(:,3), 'aa3', aa3s};
labels={labels{:},'ASB'};
longlabs={longlabs{:},'ASB - Ackerman et al.'};




% Xa. DICE with catastrophic damages separated out
index=index+1;

aa3 = defp.aa3;
aa2 = defp.aa2;
aa2_gradual = (1/(1-0.0061) - 1) / 2.5^aa3;

aa_cat_damage = 0.3; % fraction of GDP damaged in catastrophe
aa2_cat_damage = aa_cat_damage/(1-aa_cat_damage); % coefficient for the denominator

aa2_gradual_alt = (1./(1-(0.0061*calib/calib(1))) - 1) ./ repmat(2.5^aa3,1,length(calib));

cat_prob = @(Temperature) (1./(1+aa2.*abs(Temperature).^aa3) - 1./(1+aa2_gradual.*abs(Temperature).^aa3))./(-1./(1+aa2_gradual.*abs(Temperature).^aa3) + 1./(1+aa2_gradual.*abs(Temperature).^aa3 + aa2_cat_damage)); 	

cat_prob_low = @(Temperature) (1./(1+aa2*calib(2)/calib(1).*abs(Temperature).^aa3) - 1./(1+aa2_gradual_alt(:,2).*abs(Temperature).^aa3))./(-1./(1+aa2_gradual_alt(:,2).*abs(Temperature).^aa3) + 1./(1+aa2_gradual_alt(:,2).*abs(Temperature).^aa3 + aa2_cat_damage)); 	

cat_prob_high = @(Temperature) (1./(1+aa2*calib(3)/calib(1).*abs(Temperature).^aa3) - 1./(1+aa2_gradual_alt(:,3).*abs(Temperature).^aa3))./(-1./(1+aa2_gradual_alt(:,3).*abs(Temperature).^aa3) + 1./(1+aa2_gradual_alt(:,3).*abs(Temperature).^aa3 + aa2_cat_damage)); 	


% ICDF for temperature threshold for catastrophe	
testT = 0:.1:50;
testT_cat_prob=cat_prob(testT);
testT_cat_dam_multiplier = 1+ ( (testT_cat_prob>1).*(testT_cat_prob-1));

testT_cat_prob_low=cat_prob_low(testT);
testT_cat_dam_multiplier_low = 1+ ( (testT_cat_prob_low>1).*(testT_cat_prob_low-1));

testT_cat_prob_high=cat_prob_high(testT);
testT_cat_dam_multiplier = 1+ ( (testT_cat_prob_high>1).*(testT_cat_prob_high-1));


cat_thresholds = interp1(testT_cat_prob,testT,seeds(:,2),'linear');
cat_thresholds_low = interp1(testT_cat_prob_low,testT,seeds(:,2),'linear');
cat_thresholds_high = interp1(testT_cat_prob_high,testT,seeds(:,2),'linear');
		
spec_dam{index,1} = {'aa3',aa3,'aa2_gradual',aa2_gradual,'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*Temperature.^p.aa3) - 1./(1+p.aa2_gradual.*Temperature.^p.aa3))./(-1./(1+p.aa2_gradual.*Temperature.^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds,'aa_cat_width',.1,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*Temperature.^p.aa3 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)'};
spec_dam{index,2} = {'aa3',aa3,'aa2_gradual',aa2_gradual_alt(2),'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*Temperature.^p.aa3) - 1./(1+p.aa2_gradual.*Temperature.^p.aa3))./(-1./(1+p.aa2_gradual.*Temperature.^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds_low,'aa_cat_width',.1,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*Temperature.^p.aa3 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)'};
spec_dam{index,3} = {'aa3',aa3,'aa2_gradual',aa2_gradual_alt(3),'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*Temperature.^p.aa3) - 1./(1+p.aa2_gradual.*Temperature.^p.aa3))./(-1./(1+p.aa2_gradual.*Temperature.^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds_high,'aa_cat_width',.1,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*Temperature.^p.aa3 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)'};
labels={labels{:},'Xa'};
longlabs={longlabs{:},'Xa - stochastic DICE'};

% Xaa. DICE with catastrophic damages separated out and adaptation
index=index+1;

aa4 = -0.4;

spec_dam{index,1} = {'calibincomepc',BaseCalibIncomePC,'aa4',aa4, 'aa3',aa3,'aa2_gradual',aa2_gradual,'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*Temperature.^p.aa3) - 1./(1+p.aa2_gradual.*Temperature.^p.aa3))./(-1./(1+p.aa2_gradual.*Temperature.^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds,'aa_cat_width',.1,'calcdamages','1-1./(1 + (p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*Temperature.^p.aa3).*(Output./Population/p.calibincomepc).^p.aa4 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)'};
spec_dam{index,2} = {'calibincomepc',BaseCalibIncomePC,'aa4',aa4, 'aa3',aa3,'aa2_gradual',aa2_gradual_alt(2),'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*Temperature.^p.aa3) - 1./(1+p.aa2_gradual.*Temperature.^p.aa3))./(-1./(1+p.aa2_gradual.*Temperature.^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds_low,'aa_cat_width',.1,'calcdamages','1-1./(1 + (p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*Temperature.^p.aa3).*((Output./Population)/p.calibincomepc).^p.aa4 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)'};
spec_dam{index,3} = {'calibincomepc',BaseCalibIncomePC,'aa4',aa4 'aa3',aa3,'aa2_gradual',aa2_gradual_alt(3),'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*Temperature.^p.aa3) - 1./(1+p.aa2_gradual.*Temperature.^p.aa3))./(-1./(1+p.aa2_gradual.*Temperature.^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds_high,'aa_cat_width',.1,'calcdamages','1-1./(1 + (p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*Temperature.^p.aa3).*((Output./Population)/p.calibincomepc).^p.aa4 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)'};
labels={labels{:},'Xaa'};
longlabs={longlabs{:},'Xaa - stochastic DICE with adaptation'};


% Xb. DICE with catastrophic damages, eco-consumption, and capital damages separated out, and adaptation
index=index+1;

ecoshare = 0.15; % noting that at 2.5 C, 0.14% of damages (out of 0.61% gradual damages) are attributed to ecosystems/settlements, and assuming 2/3 of that represents ecosystem services share
ecoelasticity = 0.5;
dam_fcapital = 0.15; % assuming about 30% of coastal and 30% of ecosystem/settlement damages are to capital

spec_dam{index,1} = {'dam_fcapital',dam_fcapital,'ecoshare',ecoshare,'ecoelasticity',ecoelasticity,spec_dam{index-1,1}{:}};
spec_dam{index,2} = {'dam_fcapital',dam_fcapital,'ecoshare',ecoshare,'ecoelasticity',ecoelasticity,spec_dam{index-1,2}{:}};
spec_dam{index,3} = {'dam_fcapital',dam_fcapital,'ecoshare',ecoshare,'ecoelasticity',ecoelasticity,spec_dam{index-1,3}{:}};
labels={labels{:},'Xb'};
longlabs={longlabs{:},'Xb - composite'};

% Xau and Xbu. Adding uncertainty in key damage function parameters
index=index+1;

	icdfnorm = @(x,mu,v) (-sqrt(2).*erfcinv(2*x).*sqrt(v) + mu);

	N2=5000;
	clear seq;
	for j=1:100
		seq(:,1,j)=1:N2;
		for i=2:9
			seq(:,i,j)=randperm(N2);
		end
		cor = corrcoef(seq(:,:,j)).^2;
		cors(j) = sum(cor(:));
	end
	[m,mi]=min(cors);
	seeds = seq(:,:,mi)/N2 - rand(size(seq,1),size(seq,2))/N2;
	cses2 = icdfRoeBaker(seeds(:,1));

	specCS2 = {'T2xCO2',cses2};


	aa3_mean = 2;
	aa3_logvar = log(2)/2;
	aa3_mean_lognormal = log(aa3_mean) - aa3_logvar/2;
	aa3 = exp(icdfnorm(seeds(:,2),aa3_mean_lognormal,aa3_logvar));

	aa2_mean = defp.aa2;
	aa2_logvar = log(2);
	aa2_mean_lognormal = log(aa2_mean) - aa2_logvar/2;
	aa2 = exp(icdfnorm(seeds(:,7),aa2_mean_lognormal,aa2_logvar));
	aa2_gradual = (1./(1-0.0061*(aa2/defp.aa2)) - 1) ./ 2.5.^2;
	
	aa2 = aa2 .* 2.5.^(defp.aa3-aa3);
	aa2_gradual = aa2_gradual .* 2.5.^(defp.aa3-aa3);

	aa_cat_damage = 0.3 + icdftri(seeds(:,3),[-0.15 0 0.15]); 
	%aa_cat_damage = -sqrt(2).*erfcinv(2*seeds(:,3)).*0.1 + 0.3; % fraction of GDP damaged in catastrophe
	aa2_cat_damage = aa_cat_damage./(1-aa_cat_damage); % coefficient for the denominator
	
	ecoelasticity = 0.5+0.5*icdftri(seeds(:,4),[0 0.5 1]);
	ecoshare = 0.15*icdftri(seeds(:,5),[0 1 2]);;
	
	% ICDF for temperature threshold for catastrophe
	
	cat_prob = @(Temperature) (1./(1+aa2.*abs(Temperature).^aa3) - 1./(1+aa2_gradual.*abs(Temperature).^aa3))./(-1./(1+aa2_gradual.*abs(Temperature).^aa3) + 1./(1+aa2_gradual.*abs(Temperature).^aa3 + aa2_cat_damage)); 
	
	testT = 0:.1:50;
	clear testT_cat_prob;
	for i=1:length(testT)
		testT_cat_prob(:,i)=cat_prob(testT(i));
	end

	clear cat_thresholds;
	for i=1:N2
		[u,ui]=unique(testT_cat_prob(i,:));
		cat_thresholds(i) = interp1(testT_cat_prob(i,ui),testT(ui),seeds(i,6),'linear','extrap');
	end
	
	cat_thresholds=cat_thresholds(:);
	
	% fraction of damages accruing to capital
	
	dam_fcapital = 0.15 * icdftri(seeds(:,8),[0 1 2]);
	
	% adaptation
	
	aa4 = -0.4*icdftri(seeds(:,9),[0 1 2]);
	

	spec_dam{index,1} = {'aa3',aa3,'aa2_gradual',aa2_gradual,'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*abs(Temperature).^p.aa3) - 1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3))./(-1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds,'aa_cat_width',.1,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*abs(Temperature).^p.aa3 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)','T2xCO2',cses2};
	spec_dam{index,2} = spec_dam{index,1};
	spec_dam{index,3} = spec_dam{index,1};
	labels={labels{:},'Xau'};
longlabs={longlabs{:},'Xau - unc.'};

index=index+1;



	spec_dam{index,1} = {'calibincomepc',BaseCalibIncomePC,'aa4',aa4,'aa3',aa3,'aa2_gradual',aa2_gradual,'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*abs(Temperature).^p.aa3) - 1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3))./(-1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds,'aa_cat_width',.1,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*abs(Temperature).^p.aa3.*((Output./Population)/p.calibincomepc).^p.aa4 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)','T2xCO2',cses2};
	spec_dam{index,2} = spec_dam{index,1};
	spec_dam{index,3} = spec_dam{index,1};
	labels={labels{:},'Xaau'};
longlabs={longlabs{:},'Xaau - unc.'};


index = index+1;
	
	spec_dam{index,1} = {'calibincomepc',BaseCalibIncomePC,'aa4',aa4,'dam_fcapital',dam_fcapital,'ecoshare',ecoshare,'ecoelasticity',ecoelasticity,'aa3',aa3,'aa2_gradual',aa2_gradual,'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*abs(Temperature).^p.aa3) - 1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3))./(-1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3) + 1./(1+p.aa2_gradual.*Temperature.^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds,'aa_cat_width',.1,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*abs(Temperature).^p.aa3.*((Output./Population)/p.calibincomepc).^p.aa4 + (1./(1 + exp(-(Temperature-p.aa_cat_threshold)./p.aa_cat_width) )).*p.aa2_cat_damage)','T2xCO2',cses2};
	spec_dam{index,2} = spec_dam{index,1};
	spec_dam{index,3} = spec_dam{index,1};
	labels={labels{:},'Xbu'};
longlabs={longlabs{:},'Xbu - comp. unc.'};


% now calculate SCCs

clear p pCS3;
for i=1:size(spec_dam,1)
for j=1:length(etas)
for k=1:size(spec_dam,2)
for l=1:length(minoutpcs);
	p{i,j,k,l} = DICEParameters(specMCam{:},specCS{:},spec_disc{j}{:},spec_dam{i,k}{:},spec_minout{l}{:});
%	pCS3{i,j,k,l} = DICEParameters(specMCam{:},specCS{:},spec_disc{j}{:},spec_dam{i,k}{:},'T2xCO2',ones(size(cses))*3);
end
end
end
end


for t=1:10
	disp(t);
	for j=1:length(etas)
		l=1;
		k=1;
		clear a b;
		parfor i=1:size(spec_dam,1)
			disp([num2str(t) '-' num2str(i) '-' num2str(j) '-' num2str(k) '-' num2str(l)]);
			a(i) = SCC(p{i,j,k,l},[],[],t);
			b(i) = SCC(p{i,j,k,l},[],StabAbate,t);
	%		SCCs(t,i,j,k,l) = SCC(p{i,j,k,l},[],[],t);
	%		SCCStabs(t,i,j,k,l) = SCC(p{i,j,k,l},[],StabAbate,t);
		end
		SCCs(t,:,j,k,l) = a;
		SCCStabs(t,:,j,k,l) = b;
	end
end

t=2;
disp(t);
for j=1:length(etas)
	l=1;
	for k=2:size(spec_dam,2)
		clear a b;
		parfor i=1:size(spec_dam,1)
			disp([num2str(t) '-' num2str(i) '-' num2str(j) '-' num2str(k) '-' num2str(l)]);
			a(i) = SCC(p{i,j,k,l},[],[],t);
			b(i) = SCC(p{i,j,k,l},[],StabAbate,t);
	%		SCCs(t,i,j,k,l) = SCC(p{i,j,k,l},[],[],t);
	%		SCCStabs(t,i,j,k,l) = SCC(p{i,j,k,l},[],StabAbate,t);
		end
		SCCs(t,:,j,k,l) = a;
		SCCStabs(t,:,j,k,l) = b;
	end
end

t=2;
disp(t);
for j=1:length(etas)
	l=1;
	k=1;
	for l=2:length(minoutpcs)
		clear a b;
		parfor i=1:size(spec_dam,1)
			disp([num2str(t) '-' num2str(i) '-' num2str(j) '-' num2str(k) '-' num2str(l)]);
			a(i) = SCC(p{i,j,k,l},[],[],t);
			b(i) = SCC(p{i,j,k,l},[],StabAbate,t);
	%		SCCs(t,i,j,k,l) = SCC(p{i,j,k,l},[],[],t);
	%		SCCStabs(t,i,j,k,l) = SCC(p{i,j,k,l},[],StabAbate,t);
		end
		SCCs(t,:,j,k,l) = a;
		SCCStabs(t,:,j,k,l) = b;
	end
end

%save SCCs SCCs SCCsCS3 SCCStabs p spec_dam etas rhos spec_disc labels longlabs
save SCCs2 SCCs SCCStabs p spec_dam etas rhos spec_disc labels longlabs

subunc=(index-2):index;

% calculate base and ref for each
j=1;
l=1;
for k=1:size(spec_dam,2)
	clear a b
	parfor i=1:size(spec_dam,1)
		disp([num2str(i) '-' num2str(j) '-' num2str(k) '-' num2str(l)]);
		a(i) = matDICE(p{i,j,k,l}.spec{:},'calcdamages','0*Temperature');
		b(i) = matDICE(p{i,j,k,l}.spec{:});
	end
	Ref1(:,j,k,l) = a;
	Base1(:,j,k,l) = b;
end

for i=1:length(spec_dam)
disp(i)
j=1;
for k=1:size(spec_dam,2)

	Ref1(i,j,k,l).Growth = Ref1(i,j,k,l).ConsumptionPerCapita./repmat(Ref1(i,j,k,l).ConsumptionPerCapita(:,1),1,size(Ref1(i,j,k,l).ConsumptionPerCapita,2));
	Ref1(i,j,k,l).EffectiveGrowth = Ref1(i,j,k,l).EffectiveConsumptionPerCapita./repmat(Ref1(i,j,k,l).EffectiveConsumptionPerCapita(:,1),1,size(Ref1(i,j,k,l).EffectiveConsumptionPerCapita,2));

	Base1(i,j,k,l).Growth = Base1(i,j,k,l).ConsumptionPerCapita./repmat(Base1(i,j,k,l).ConsumptionPerCapita(:,1),1,size(Base1(i,j,k,l).ConsumptionPerCapita,2));
	Base1(i,j,k,l).EffectiveGrowth = Base1(i,j,k,l).EffectiveConsumptionPerCapita./repmat(Base1(i,j,k,l).EffectiveConsumptionPerCapita(:,1),1,size(Base1(i,j,k,l).EffectiveConsumptionPerCapita,2));
end
end

clear DiscRateMedian_Ref1 DiscRateMean_Ref1 DiscRateMedian_Base1 DiscRateMean_Base1
for i=1:length(spec_dam)
disp(i)
for k=1:size(spec_dam,2)
for j=1:length(etas)
l=1;
	DiscFact_Ref = repmat(p{i,j,k,l}.rr,size(Ref1(i,1,k).Growth,1),1).*(Ref1(i,1,k).Growth.^-p{i,j,k}.elasmu);
	EffDiscFact_Ref = repmat(p{i,j,k,l}.rr,size(Ref1(i,1,k).Growth,1),1).*(Ref1(i,1,k).EffectiveGrowth.^-p{i,j,k,l}.elasmu);


	DiscRateMedian_Ref1(i,:,j,k) = quantile((repmat(DiscFact_Ref(:,2),1,size(DiscFact_Ref,2))./DiscFact_Ref).^(1./(repmat(Ref1(i,1,k).times,size(DiscFact_Ref,1),1)-Ref1(i,1,k).times(2)+eps)),.5)-1;
	DiscRateMean_Ref1(i,:,j,k) = mean((repmat(DiscFact_Ref(:,2),1,size(DiscFact_Ref,2))./DiscFact_Ref).^(1./(repmat(Ref1(i,1,k).times,size(DiscFact_Ref,1),1)-Ref1(i,1,k).times(2)+eps)))-1;

	DiscFact_Base = repmat(p{i,j,k,l}.rr,size(Base1(i,1,k).Growth,1),1).*(Base1(i,1,k).Growth.^-p{i,j,k,l}.elasmu);
	EffDiscFact_Base = repmat(p{i,j,k,l}.rr,size(Base1(i,1,k).Growth,1),1).*(Base1(i,1,k).EffectiveGrowth.^-p{i,j,k,l}.elasmu);


	DiscRateMedian_Base1(i,:,j,k) = quantile((repmat(DiscFact_Base(:,2),1,size(DiscFact_Base,2))./DiscFact_Base).^(1./(repmat(Base1(i,1,k).times,size(DiscFact_Base,1),1)-Base1(i,1,k).times(2)+eps)),.5)-1;
	DiscRateMean_Base1(i,:,j,k) = mean((repmat(DiscFact_Base(:,2),1,size(DiscFact_Base,2))./DiscFact_Base).^(1./(repmat(Base1(i,1,k).times,size(DiscFact_Base,1),1)-Base1(i,1,k).times(2)+eps)))-1;


end
end
end


%%%%% START FIGURES
%%%%%

% Table: Discount Rates

m=1;
headA = {};
for i=1:length(etas)
	headA = {headA{:},['rho = ' num2str(rhos(i)*100) ' pct / eta = ' num2str(etas(i)) ]};
end

headB = {'median to 2115','expected to 2115','median to 2215','expected to 2215'};
headerrow=[];
for i=1:length(headA)
for j=1:length(headB)
headerrow=[headerrow headA{i} ' - ' headB{j} ','];
end
end
labelrow = 'Specs: Ref';
clear DiscTable;
ind = 1;
for k=1:length(etas)
	DiscTable(ind,(k-1)*4+1) = DiscRateMedian_Ref1(1,12,k,1);
	DiscTable(ind,(k-1)*4+2) = DiscRateMean_Ref1(1,12,k,1);
	DiscTable(ind,(k-1)*4+3) = DiscRateMedian_Ref1(1,22,k,1);
	DiscTable(ind,(k-1)*4+4) = DiscRateMean_Ref1(1,22,k,1);
end
ind = ind + 1;
for j=1:length(spec_dam)
	labelrow=[labelrow ',' labels{j}];
	for k=1:length(etas)
		DiscTable(ind,(k-1)*4+1) = DiscRateMedian_Base1(j,12,k,1);
		DiscTable(ind,(k-1)*4+2) = DiscRateMean_Base1(j,12,k,1);
		DiscTable(ind,(k-1)*4+3) = DiscRateMedian_Base1(j,22,k,1);
		DiscTable(ind,(k-1)*4+4) = DiscRateMean_Base1(j,22,k,1);
	end
	ind = ind+1;
end

filename='DiscTable.csv'; towrite=DiscTable*100;
fid=fopen(filename,'w');
fprintf(fid,[labelrow '\n\n' headerrow '\n']);
fclose(fid);
dlmwrite(filename,towrite,'delimiter',',','precision','%.2f','-append');

%%

% change in deterministic SCC with discount rate


%for i=1:length(spec_dam)
%for k=1:length(etas)-1
%	ratdetSCC(i,:,k) = [SCCs(2,i,k+1,1,1).SCC_consumptiondenominated]'./[SCCs(2,i,1,1,1).SCC_consumptiondenominated]';
%	[ratdetsort(i,:,k),ratdetsorti(i,:,k)] = sort(ratdetSCC(i,:,k));
%end
%end	
%u=ratdetsorti(:,500,1);
%for i=1:length(u)
%	ratdetmid(i,:) = squeeze(ratdetsort(i,u(i),:));
%end

%%


colors = 'bgrcmk';
lines = {'-','--','-.'};

% Figure: Reference Scenarios

figure;
clf;
clear ha;
ha(1)=subplot(2,2,1)
plot(p{1,1,1,1}.t,mean([RefUnc(1).Emissions])/10*(44/12),colors(1),'LineWidth',2); hold on
plot(p{1,1,1,1}.t,mean([RefStabUnc(1).Emissions])/10*(44/12),colors(2),'LineWidth',2);
plot(p{1,1,1,1}.t,quantile([RefUnc(1).Emissions],.5)/10*(44/12),colors(1));
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).Emissions],.5)/10*(44/12),colors(2));
plot(p{1,1,1,1}.t,quantile([RefUnc(1).Emissions],[.05 .95])/10*(44/12),[colors(1) '--']);
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).Emissions],[.05 .95])/10*(44/12),[colors(2) '--']);
xlim([2000 2300]);
xlabel('Year');ylabel('Emissions (Gt CO_2/year)')


ha(2)=subplot(2,2,2)
plot(p{1,1,1,1}.t,mean([RefUnc(1).ppmCO2]),colors(1),'LineWidth',2); hold on
plot(p{1,1,1,1}.t,mean([RefStabUnc(1).ppmCO2]),colors(2),'LineWidth',2);
plot(p{1,1,1,1}.t,quantile([RefUnc(1).ppmCO2],.5),colors(1));
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).ppmCO2],.5),colors(2));
plot(p{1,1,1,1}.t,quantile([RefUnc(1).ppmCO2],[.05 .95]),[colors(1) '--']);
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).ppmCO2],[.05 .95]),[colors(2) '--']);
plot(p{1,1,1,1}.t,mean(278*2.^[RefUnc(1).Forcing/defp.FCO22x]),[colors(1) '-.']); hold on
plot(p{1,1,1,1}.t,mean(278*2.^[RefStabUnc(1).Forcing/defp.FCO22x]),[colors(2) '-.']);
%legend('Ref','Base','Optim','Panic','\Delta Panic','Location','NorthWest');
xlim([2000 2300]);
xlabel('Year');ylabel('CO_2 or CO_2e conc. (ppm)')

ha(3)=subplot(2,2,3)
plot(p{1,1,1,1}.t,mean([RefUnc(1).Tatm]),colors(1),'LineWidth',2); hold on
plot(p{1,1,1,1}.t,mean([RefStabUnc(1).Tatm]),colors(2),'LineWidth',2);
plot(p{1,1,1,1}.t,quantile([RefUnc(1).Tatm],.5),colors(1));
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).Tatm],.5),colors(2));
plot(p{1,1,1,1}.t,quantile([RefUnc(1).Tatm],[.05 .95]),[colors(1) '--']);
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).Tatm],[.05 .95]),[colors(2) '--']);
xlim([2000 2300]);
xlabel('Year');ylabel('T_{atm} (^oC)')

ha(4)=subplot(2,2,4)
plot(p{1,1,1,1}.t,mean([RefUnc(1).Output]),colors(1),'LineWidth',2); hold on
plot(p{1,1,1,1}.t,mean([RefStabUnc(1).Output]),colors(2),'LineWidth',2);
plot(p{1,1,1,1}.t,quantile([RefUnc(1).Output],.5),colors(1));
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).Output],.5),colors(2));
plot(p{1,1,1,1}.t,quantile([RefUnc(1).Output],[.05 .95]),[colors(1) '--']);
plot(p{1,1,1,1}.t,quantile([RefStabUnc(1).Output],[.05 .95]),[colors(2) '--']);
hl = legend('Ref.','Stab.','Location','SouthEast');
xlim([2000 2300]);
xlabel('Year');ylabel('GDP (Trillion USD)')

[bh,th]=label(ha,'ul',12,[],0,1,1,1.5,1.5); delete(bh);
longticks
pdfwrite(gcf,'Scenarios');


% Figure: expected damages

% come up with expected damages

testT = [0:.1:10]';

clear dam damq dam_effs dam_effsq damA dam_effsA damqA dam_effsqA;
j=1; l=1;
for i=1:size(spec_dam,1)
	disp(i);
	for k=1:size(spec_dam,2)
		ecosh = p{i,j,k,l}.ecoshare;
		ecoel = p{i,j,k,l}.ecoelasticity;
		d1 = p{i,j,k,l}.calcdamages(p{i,j,k,l}.Tatm0,p{i,j,k,l}.t(1),p{i,j,k,l}.Trate0,p{i,j,k,l}.q0,p{i,j,k,l}.Tatm0 - p{i,j,k,l}.Trate0 * [-100:10:0],p{i,j,k,l}.pop0);
		
		for n=1:length(testT)
			%q = p{i,j,k,l}.q0;
			%rate = p{i,j,k,l}.Trate0;
			q = BaseCalibGDP;
			pop = BaseCalibPop;
			rate = BaseCalibRate;
			%d = p{i,j,k,l}.calcdamages(testT(j),p{i,j,k,l}.t(1),rate*10,q,repmat(testT(j),1,11) + repmat(rate * [-100:10:0],length(testT(j)),1));
			d = p{i,j,k,l}.calcdamages(testT(n),BaseCalibTime,rate*10,q,repmat(testT(n),1,11) + repmat(BaseCalibRate * [-100:10:0],length(testT(n)),1),pop);
			dam(n,i,j,k,l) = mean(d);
			damq(n,i,j,k,l,:) = reshape(quantile(d,[.05 .5 .95 .01 .99],1),1,1,1,[]);
			dam_effs(n,i,j,k,l)= mean ( (q - ( ( ( (1-ecosh).*(q * (1-d)).^(1-1./ecoel) + ecosh .* (p{i,j,k,l}.q0 * (1-d)).^(1-1./ecoel) - ecosh .* (p{i,j,k,l}.q0 * (1-d1)).^(1-1./ecoel(:,1)) ) ./ (1-ecosh)) .^ (ecoel./(ecoel-1)) ) ) / q);
			dam_effsq(n,i,j,k,l,:) = reshape(quantile( (q - ( ( ( (1-ecosh).*(q * (1-d)).^(1-1./ecoel) + ecosh .* (p{i,j,k,l}.q0 * (1-d)).^(1-1./ecoel) - ecosh .* (p{i,j,k,l}.q0 * (1-d1)).^(1-1./ecoel(:,1)) ) ./ (1-ecosh)) .^ (ecoel./(ecoel-1)) ) ) / q,[.05 .5 .95 .01 .99],1),1,1,1,[]);

			dA = p{i,j,k,l}.calcdamages(testT(n),BaseCalibTime,rate*10*2,q*2,repmat(testT(n),1,11) + repmat(BaseCalibRate*2 * [-100:10:0],length(testT(n)),1),pop);
			damA(n,i,j,k,l) = mean(dA);
			damqA(n,i,j,k,l,:) = reshape(quantile(d,[.05 .5 .95 .01 .99],1),1,1,1,[]);
			dam_effsA(n,i,j,k,l)= mean ( (q - ( ( ( (1-ecosh).*(q * (1-d)).^(1-1./ecoel) + ecosh .* (p{i,j,k,l}.q0 * (1-d)).^(1-1./ecoel) - ecosh .* (p{i,j,k,l}.q0 * (1-d1)).^(1-1./ecoel(:,1)) ) ./ (1-ecosh)) .^ (ecoel./(ecoel-1)) ) ) / q);
			dam_effsqA(n,i,j,k,l,:) = reshape(quantile( (q - ( ( ( (1-ecosh).*(q * (1-d)).^(1-1./ecoel) + ecosh .* (p{i,j,k,l}.q0 * (1-d)).^(1-1./ecoel) - ecosh .* (p{i,j,k,l}.q0 * (1-d1)).^(1-1./ecoel(:,1)) ) ./ (1-ecosh)) .^ (ecoel./(ecoel-1)) ) ) / q,[.05 .5 .95 .01 .99],1),1,1,1,[]);


		end
	end
end


subdistinct = [1 2 3 7 8 9 13];
subdistinctA = [3 5 6 11 14];

maxdam95 = squeeze(max(damq(:,subdistinct(1:end-1),j,:,l,3),[],2));
mindam5 = squeeze(min(damq(:,subdistinct(1:end-1),j,:,l,1),[],2));
maxdam_eff95 = squeeze(max(dam_effsq(:,subdistinct(1:end-1),j,:,l,3),[],2));
mindam_eff5 = squeeze(min(dam_effsq(:,subdistinct(1:end-1),j,:,l,1),[],2));


j=1; k=1; l=1;
figure;
clear ha;
ha(1)=subplot(2,1,1)
hpb=plot(testT,[mindam5(:,k) maxdam95(:,k)]*100,'--','Color',[.6 .6 .6],'LineWidth',2); hold on;
hp = plot(testT,dam(:,subdistinct,j,k,l)*100);
hl = legend(hp,labels{subdistinct},'Location','NorthWest')
for i=1:length(hp)
	set(hp(i),'Color',colors(1+mod(i-1,length(colors))),'LineStyle',lines{ceil((i)/length(colors))});
end
set(hp(1),'linewidth',2);

ylim([0 30]);
ht=text(7.5,4,{'dT/dt = 0.37^oC/dec','C = $12.6K/person/yr'});
%title('Damages - Initial Output');
ylabel({'Expected Damages (%)'}); xlabel('Warming (^oC)');
%ps=get(gca,'position'); ps(3)=0.55; set(gca,'position',ps);
%pdfwrite('damages');

%figure;
ha(2)=subplot(2,1,2)
hp = plot(testT,dam(:,subdistinct(1),j,k,l)*100); hold on;
hpA = plot(testT,damA(:,subdistinctA,j,k,l)*100);
clear leg;
leg{1}='D';
for i=1:length(subdistinctA)
	leg{end+1} = [labels{subdistinctA(i)} char(39)];
end
hl = legend([hp ; hpA],leg,'Location','NorthWest')
for i=1:length(hp)
	set(hp(i),'Color',colors(1+mod(i-1,length(colors))),'LineStyle',lines{ceil((i)/length(colors))});
end
colorsA = 'brgcmk';
for i=1:length(hpA)
	set(hpA(i),'Color',colorsA(1+mod(i-1+length(hp),length(colors))),'LineStyle',lines{ceil((i+length(hp))/length(colors))});
end
set(hp(1),'linewidth',2);
hold on
%hp(end+1) = plot(testT,dam5xGDP(:,3,1)*100,':','Color',get(hp(3),'Color'))
%hp(end+1) = plot(testT,dam5xGDP(:,5,1)*100,':','Color',get(hp(5),'Color'))

ylim([0 30]);
ht=text(7.5,4,{'dT/dt = 0.74^oC/dec','C = $25.2K/person/yr'});
ylabel({'Expected Damages (%)'}); xlabel('Warming (^oC)');
%ps=get(gca,'position'); ps(3)=0.55; set(gca,'position',ps);

[bh,th]=label(ha,'ll',12,[],0,1,1,1.5,1.5); %delete(bh);
ylabel(ha(1),''); movev(ha(2),.05);
pdfwrite('damages');

clf;
subplot(2,1,1)
hp0=plot(testT,dam(:,1,j,k,l)*100,'b--','LineWidth',2);
hold on;
hpa=plot(testT,[squeeze(damq(:,end,j,k,l,1)) squeeze(damq(:,end,j,k,l,3))]*100,'k--');
hpb=plot(testT,[squeeze(damq(:,end,j,k,l,4)) squeeze(damq(:,end,j,k,l,5))]*100,'k:');
hpc=plot(testT,[squeeze(damq(:,end,j,k,l,2))]*100,'k-');
hpc=plot(testT,dam(:,end,j,k,l)*100,'k-','linewidth',2);
ylim([0 100]);
title('Damages Distribution - \bf{Xau}, \bf{Xaau}, \bf{Xbu}');ylabel('Damages (%)'); xlabel('Warming (^oC)');
pdfwrite('damages_Xbu');

% Figure: GDP in 2115

clf;
subplot(2,1,1)
for i=1:size(p,1)
	u = [Base1(i,j,k,l).ConsumptionPerCapita]; uref=[Ref1(i,j,k,l).ConsumptionPerCapita];
	u2 = [Base1(i,j,k,l).EquivalentMaterialConsumptionPerCapita]; uref2=[Ref1(i,j,k,l).EquivalentMaterialConsumptionPerCapita];
	u=u(:,12)./u(:,2); uref=uref(:,12)./uref(:,2);
	u2=u2(:,12)./u2(:,2); uref2=uref2(:,12)./uref2(:,2);
	
%	plot([i i i]-.1,quantile(uref,[.05 .5 .95]),'k.-'); hold on
	plot([i i i],quantile(u,[.05 .5 .95]),'b.-'); hold on
	if mean(abs(u-u2))>.05
%		plot([i i i]+.2,quantile(uref2,[.05 .5 .95]),'k.--'); hold on
		plot([i i i]+.3,quantile(u2,[.05 .5 .95]),'b.--'); hold on
	end
%	plot([i i i]-.1,mean(uref),'kx'); hold on
	plot(i,mean(u),'bs'); hold on
	if mean(abs(u-u2))>.05
%		plot([i i i]+.2,mean(uref2),'kx'); hold on
		plot(i+.3,mean(u2),'bd'); hold on
	end
end
%plot(0,mean(uref),'ks');
plot([-10 100],[mean(uref) mean(uref)],'k:');
xlim([0 size(p,1)+1]); ylim([3 8]);
set(gca,'XTickLabel',{labels{:}},'XTick',1:length(labels));
ylabel({'2115 Consumption Per Capita','(2015 = 1)'});
longticks(gca,2)
pdfwrite('ConsPerCap');


% Figure: GDP path
j=1;k=1;l=1;
for i=1:length(spec_dam)
	clf;
	u = [Base1(i,j,k,l).ConsumptionPerCapita]; uref=[Ref1(i,j,k,l).ConsumptionPerCapita];
	u2 = [Base1(i,j,k,l).EquivalentMaterialConsumptionPerCapita]; uref2=[Ref1(i,j,k,l).EquivalentMaterialConsumptionPerCapita];
	u=u./repmat(u(:,2),1,size(u,2)); uref=uref./repmat(uref(:,2),1,size(u,2));
	u2=u2./repmat(u2(:,2),1,size(u,2)); uref2=uref2./repmat(uref2(:,2),1,size(u,2));
	
	subplot(2,1,1)
	plot(Base1(i,j,k,l).times,mean(uref),'k-'); hold on
	plot(Base1(i,j,k,l).times,quantile(u,[.05 .95]),'b--'); hold on
	plot(Base1(i,j,k,l).times,quantile(u,.5),'b-'); hold on
	plot(Base1(i,j,k,l).times,mean(u),'b-','LineWidth',2); hold on
	ylabel('Consumption Per Capita (2015 = 1)');
	xlim([2000 2300]);
	title(['\bf{' labels{i} '}']);
	
	
%	
%	if mean(abs(u(:,15)-u2(:,15)))>.05
%		subplot(2,1,2)
%		plot(Base1(i).times,mean(uref2),'k-'); hold on
%		plot(Base1(i).times,quantile(u2,[.05 .5 .95]),'b--'); hold on
%		plot(Base1(i).times,quantile(u2,.5),'b-'); hold on
%		plot(Base1(i).times,mean(u2),'b-','LineWidth',2); hold on
%		ylabel('Effective Consumption Per Capita (2015 = 1)');
%	end
	pdfwrite(['ConsPerCapPath_' labels{i}]);
end

% Figure: SCCs

clf;
k=1; l=1;
clear ha sub; sub{1} = [1:9]; sub{2} = [1 10:15];
for m=1:2
	ha(m) = subplot(2,1,m)
	for i=1:length(sub{m})
	for j=1:length(etas)
		plot([i i i]+(j-2.5)*.2,quantile(SCCs(2,sub{m}(i),j,k,l).SCC_materialconsumptionequiv,[.05 .5 .95]),[colors(j) '.-']); hold on
	end
	end
	for i=1:length(sub{m})
	for j=1:length(etas)
		plot(i+(j-2.5)*.2,SCCs(2,sub{m}(i),j,k,l).ESCC_materialconsumptionequiv,[colors(j) 's']); hold on
	end
	end
	set(gca,'XTickLabel',labels(sub{m}),'XTick',1:length(sub{m}));
	ylabel('2015 SCC ($/ton CO_2)');
	leg = {};
end
for j=1:length(etas)
	leg = {leg{:},['\eta=' num2str(etas(j)) ', \rho=' sprintf('%0.2g',(rhos(j)*100)) '%' ]};
end
hl=legend(leg{:},'Location','East'); set(ha,'xlim',[.5 10],'YScale','log','ylim',[0.5 2000],'YTick',10.^[-1:4],'YMinorTick','on');
ylim(ha(1),[5 2000]);
title(ha(1),'SCCs - Reference Scenario');
movev(ha(2),.07);
pdfwrite('SCCs');

clf;
k=1; l=1;
clear ha;
for m=1:2
	ha(m) = subplot(2,1,m)
	for i=1:length(sub{m})
	for j=[1 3]
		plot([i i i]+(j-2.5)*.2-.1,quantile(SCCs(2,sub{m}(i),j,k,l).SCC_materialconsumptionequiv,[.05 .5 .95]),[colors(j) '.-']); hold on
		plot([i i i]+(j-2.5)*.2+.1,quantile(SCCStabs(2,sub{m}(i),j,k,l).SCC_materialconsumptionequiv,[.05 .5 .95]),[colors(j) '.--']); hold on
	end
	end
	for i=1:length(sub{m})
	for j=[1 3]
		plot(i+(j-2.5)*.2-.1,SCCs(2,sub{m}(i),j,k,l).ESCC_materialconsumptionequiv,[colors(j) 's']); hold on
		plot(i+(j-2.5)*.2+.1,SCCStabs(2,sub{m}(i),j,k,l).ESCC_materialconsumptionequiv,[colors(j) 'd']); hold on
	end
	end
	set(gca,'XTickLabel',labels(sub{m}),'XTick',1:length(sub{m}));
	ylabel('2015 SCC ($/ton CO_2)');
end
leg = {};
for j=[1 3]
	disclang = ['\eta=' num2str(etas(j)) ', \rho=' sprintf('%0.2g',(rhos(j)*100)) '%' ];
	leg = {leg{:},['Ref: ' disclang ],['Stab: ' disclang ]};
end
hl=legend(ha(2),leg{:},'Location','East');
set(hl,'fontsize',8);
  set(ha,'xlim',[.5 10],'YScale','log','ylim',[0.5 2000],'YTick',10.^[-1:4],'YMinorTick','on');
ylim(ha(1),[5 2000]);
htitle=title(ha(1),'SCCs - Reference Scenario vs. Stabilization Scenario');
movev(ha(2),.07);
pdfwrite('SCCs_RefVsStab');


% Table: SCCs - Reference

m=1; l=1;
headA = {};
for i=1:length(etas)
	headA = {headA{:},['rho = ' num2str(rhos(i)*100) ' pct / eta = ' num2str(etas(i)) ]};
end

headB = {'expected','median','5th','95th'};
headerrow=[];
for i=1:length(headA)
for j=1:length(headB)
headerrow=[headerrow headA{i} ' - ' headB{j} ','];
end
end
labelrow = 'Specs: ';
clear SCCTable;
ind = 1;
for j=1:length(spec_dam)
	labelrow=[labelrow ',' labels{j}];
	for k=1:length(etas)
		SCCTable(ind,(k-1)*4+1) = [SCCs(2,j,k,m,l).ESCC_materialconsumptionequiv];
		SCCTable(ind,(k-1)*4+[2:4]) = quantile([SCCs(2,j,k,m,l).SCC_materialconsumptionequiv],[.5 0.05 0.95]);
	end
	ind = ind+1;
	if abs(SCCs(2,j,1,m).ESCC_materialconsumptionequiv-SCCs(2,j,1,m,l).ESCC_consumptiondenominated)>.05
		labelrow=[labelrow ',' labels{j} '*'];
		for k=1:length(etas)
			SCCTable(ind,(k-1)*4+1) = [SCCs(2,j,k,m,l).ESCC_consumptiondenominated];
			SCCTable(ind,(k-1)*4+[2:4]) = quantile([SCCs(2,j,k,m,l).SCC_consumptiondenominated],[.5 0.05 0.95]);
		end
		ind = ind+1;
	end

end

filename='SCCTable.csv'; towrite=SCCTable;
fid=fopen(filename,'w');
fprintf(fid,[labelrow '\n\n' headerrow '\n']);
fclose(fid);
dlmwrite(filename,towrite,'delimiter',',','precision','%.2f','-append');


% Table: SCCs - Stabilization

m=1; l=1;
headA = {};
for i=1:length(etas)
	headA = {headA{:},['rho = ' num2str(rhos(i)*100) ' pct / eta = ' num2str(etas(i)) ]};
end

headB = {'expected','median','5th','95th'};
headerrow=[];
for i=1:length(headA)
for j=1:length(headB)
headerrow=[headerrow headA{i} ' - ' headB{j} ','];
end
end
labelrow = 'Specs: ';
clear SCCTable;
ind = 1;
for j=1:length(spec_dam)
	labelrow=[labelrow ',' labels{j}];
	for k=1:length(etas)
		SCCTable(ind,(k-1)*4+1) = [SCCStabs(2,j,k,m,l).ESCC_materialconsumptionequiv];
		SCCTable(ind,(k-1)*4+[2:4]) = quantile([SCCStabs(2,j,k,m,l).SCC_materialconsumptionequiv],[.5 0.05 0.95]);
	end
	ind = ind+1;
	if abs(SCCs(2,j,1,m).ESCC_materialconsumptionequiv-SCCs(2,j,1,m,l).ESCC_consumptiondenominated)>.05
		labelrow=[labelrow ',' labels{j} '*'];
		for k=1:length(etas)
			SCCTable(ind,(k-1)*4+1) = [SCCStabs(2,j,k,m,l).ESCC_consumptiondenominated];
			SCCTable(ind,(k-1)*4+[2:4]) = quantile([SCCStabs(2,j,k,m,l).SCC_consumptiondenominated],[.5 0.05 0.95]);
		end
		ind = ind+1;
	end

end

filename='SCCTable_Stabilization.csv'; towrite=SCCTable;
fid=fopen(filename,'w');
fprintf(fid,[labelrow '\n\n' headerrow '\n']);
fclose(fid);
dlmwrite(filename,towrite,'delimiter',',','precision','%.2f','-append');

% Figure: varying calibration parameters

clf;
subplot(2,1,1)
ind=0;
xtickl={'0'};
clear hp d;
l=1;
for j=[1]
	ind = ind + 1;
	for i=1:((subunc(1)-1))
		for k=1:3
			d(i,k) = mean(p{i,j,k,l}.calcdamages(2.5,BaseCalibTime,BaseCalibRate*10,BaseCalibGDP,2.5 + BaseCalibRate * [-100:10:0],BaseCalibPop));
		end
		hp(i)=plot([0 d(i,[2 1 3])]*100,[0 SCCs(2,i,j,[2 1 3],l).ESCC_materialconsumptionequiv]/SCCs(2,i,j,1,l).ESCC_materialconsumptionequiv,'s-','Color',colors(1+mod(i-1,length(colors))),'LineStyle',lines{ceil((i)/length(colors))}); hold on
	end
	ylabel('2015 SCC (Central Calib = 1)');
	title(['\eta = ' num2str(etas(j)) ' / \rho = ' sprintf('%0.2g',rhos(j)*100) '%']);
%	set(gca,'YScale','log','XScale','log')
	xlim([0 calib(3)*105]); set(gca,'xtick',[0 1.8 3.6]);
end
xlabel('Expected damage at 2.5^oC (%)');
hl = legend(hp,labels{1:end-1},'Location','EastOutside')
ps=get(gca,'position'); ps(3)=0.55; set(gca,'position',ps);
pdfwrite('SCCs_calibration');

% Figure: SCC paths


clear ESCC ESCCeq
m=1; l=1;
for i=1:length(spec_dam)
for k=1:length(etas)
	ESCC(i,:,k,m) = [SCCs(:,i,k,m,l).ESCC_consumptiondenominated];
	ESCCeq(i,:,k,m) = [SCCs(:,i,k,m,l).ESCC_materialconsumptionequiv];
end
end

for k=1:length(etas)
	subs={[1:4],[5:8],[9:12],[13:14]};
	clf
	for i=1:length(subs)
		subplot(length(subs),1,i)
		plot(p{1,1,1,1}.t(1:size(ESCC,2)),ESCCeq(subs{i},:,k)); hold on
		plot(p{1,1,1,1}.t(1:size(ESCC,2)),ESCC(subs{i},:,k),':');
		ylabel('SCC ($/ton CO_2)');
		 legend(labels{subs{i}});
		if i==1
			title(['\eta = ' num2str(etas(k)) ', \rho = ' sprintf('%0.2g',(rhos(k)*100)) ]);
		end		
	end
	pdfwrite(['SCCs_path' num2str(k)]);
end

clear ESCC ESCCeq
m=1; l=1;
for i=1:length(spec_dam)
for k=1:length(etas)
	ESCC(i,:,k,m,l) = [SCCStabs(:,i,k,m,l).ESCC_consumptiondenominated];
	ESCCeq(i,:,k,m,l) = [SCCStabs(:,i,k,m,l).ESCC_materialconsumptionequiv];
end
end

for k=1:length(etas)
	subs={[1:4],[5:8],[9:12],[13:14]};
	clf
	for i=1:length(subs)
		subplot(length(subs),1,i)
		plot(p{1}.t(1:size(ESCC,2)),ESCCeq(subs{i},:,k)); hold on
		plot(p{1}.t(1:size(ESCC,2)),ESCC(subs{i},:,k),':');
		ylabel('SCC ($/ton CO_2)');
		 legend(labels{subs{i}});
		if i==1
			title(['Stabilization Scenario - \eta = ' sprintf('%0.2g',(etas(k)*100)) '% / \rho = ' num2str(rhos(k)) ]);
		end		
	end
	pdfwrite(['SCCs_Stab_path' num2str(k)]);
end


% Figure: SCC path - D and Xc at two discount rates, stabilization vs. reference

clear ESCC ESCCeq
m=1; l=1;
for i=1:length(spec_dam)
for k=1:length(etas)
	ESCC(i,:,k,m,l) = [SCCs(:,i,k,m,l).ESCC_consumptiondenominated];
	ESCCeq(i,:,k,m,l) = [SCCs(:,i,k,m,l).ESCC_materialconsumptionequiv];
	ESCCStab(i,:,k,m,l) = [SCCStabs(:,i,k,m,l).ESCC_consumptiondenominated];
	ESCCeqStab(i,:,k,m,l) = [SCCStabs(:,i,k,m,l).ESCC_materialconsumptionequiv];
end
end

clf;
ind=0;
for k=[1 3]
	ind=ind+1;
	subplot(3,1,ind)
	sub=[1 14];
	plot(p{1}.t(1:size(ESCC,2)),ESCCeq(sub,:,k),'linew',2); hold on
	plot(p{1}.t(1:size(ESCC,2)),ESCCeqStab(sub,:,k),'--','linew',2); hold on
	ylabel('SCC ($/ton CO_2)');
	title(['\eta=' num2str(etas(k)) ', \rho=' sprintf('%0.2g',(rhos(k)*100)) '%' ]);
	
end
	legend('D - Ref','Xbu - Ref','D - Stab','Xbu - Stab','Location','Northwest');

pdfwrite(['SCCs_RefVsStab_path']);

%%%%%%%%%%%
%%%%%%%%%%%

% Figure: ternary diagram fcap and futil

defp=DICEParameters;

fcap = 0:.05:1;
futil = 0:.05:1;
[FUTIL,FCAP]=meshgrid(fcap,futil);
for i=1:length(fcap)
disp(i);
for j=1:length(futil)

	if (fcap(i)+futil(j)>1)
		SCCcaputilA(i,j) = NaN;
		SCCcaputilB(i,j) = NaN;
	else
		a = SCC(DICEParameters(specMCam{:},spec_disc{1}{:},'dam_fcapital',fcap(i),'dam_futility',futil(j)));
		b = SCC(DICEParameters(specMCam{:},spec_disc{1}{:},'dam_fcapital',fcap(i),'dam_futility',futil(j),'aa3',4,'aa2',defp.aa2*2.5^(2-4)));		SCCcaputilA(i,j) = a.ESCC_materialconsumptionequiv;
		SCCcaputilB(i,j) = b.ESCC_materialconsumptionequiv;		
	end

end
end

clf;
sub=find(~isnan(SCCcaputilA));
terncontourf(FCAP(sub),FUTIL(sub),SCCcaputilA(sub));
ht = ternlabel('% capital','% utility','% output'); 
%title('DICE Damage Function, 3% discount rate');
hc = colorbar('SouthOutside'); xl=xlabel(hc,'$/ton CO_2');
set(ht,'fontsize',12)
set(hc,'Position',[.25 .08 .5 .0595])
movev(ht(1),-.01);
moveh(ht(2),.04);
moveh(ht(3),-.04);
ht1(1) = text(-.04,.02,'Output');
ht1(2) = text(1.06,.02,'Capital');
ht1(3) = text(0.5,.9,'Utility');
set(ht1,'fontsize',12)
set(ht1(1),'HorizontalAlignment','Right');
set(ht1(2),'HorizontalAlignment','Left');
set(ht1(3),'HorizontalAlignment','Center');
pdfwrite('SCC_capital_utility_quadratic');

%%%%%%%
%%%%%%%

% Figure: Mappings

defp=DICEParameters;
testT = 0:.01:30;
straightmapping = defp.aa2 * testT.^2;
rationalmapping = 1 - 1./(1+straightmapping);
expmapping = 1 - exp(-straightmapping);

clear ha;
clf;
subplot(2,1,1);
ha(1)=plot(testT,straightmapping*100,'k'); hold on;
ha(2)=plot(testT,rationalmapping*100,'b--');
ha(3)=plot(testT,expmapping*100,'r-.');
ylim([0 150]);
legend('Direct','Rational (D)','Exponential (We)','Location','Northwest');
xlabel('Warming (^oC)'); ylabel('Damages (% of Output)');
set(ha,'linew',2);
longticks(gca,2);
pdfwrite('mappings');