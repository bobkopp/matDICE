function p=DICEParameters(varargin)

	% p=DICEParameters([parameter],[value],...)
	%
	% Example:
	%     p = DICEParameters('elasmu',0,'prstp',0.03);
	%
	% Default values largely match DICE 2010
	%
	% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Mar 4 19:06:31 EST 2013
	% Copyright (C) 2012 by Robert E. Kopp; distributed under GNU GPL v3

	p.spec = varargin;

	p.t=2005:10:2305;
	
	% Preferences
	
	p.elasmu=1.5; %Elasticity of marginal utility of consumption
	p.prstp=0.015; %Initial rate of social time preference per year
	
	% Population and technology
	p.pop0=6411; % 2005 world population millions
	p.popasym=8700; % Asymptotic population
	%p.gpop0=.115; % Growth rate of population per decade
	p.popadjrate=0.59; % Population adjustment rate per decade
	p.ga0=.160; % Initial growth rate of productivity per decade
	p.dela=.0094; % Decline rate of productivity growth per decade
	p.deldela=0.0019; % Rate of decline of decline rate in productivity growth
	p.dk=.100; % Depreciation rate on capital per year
	p.gama=.300; % Capital elasticity in production function
	p.q0=57.09; % 2005 world gross output trill 2005 US international dollars % matches world bank figures
	p.k0=97.3/55.3 * p.q0; %2005 value capital trill 2005 US dollars
	%p.a0=.02722; % Initial level of total factor productivity
	p.a0 = p.q0/(p.k0^p.gama * p.pop0^(1-p.gama));
	p.minoutpc = 730; % dollars per person per year
	p.mincons = 730; % dollars per person per year
	p.popfeedback = 0; % population feedback on or off (do we react to output less than the minimum by raising output or cutting population?)
	p.allowrecoveryfromcrash = 0; % are crashes to minimum output forever?
	
	% Emissions
	
	p.CO2_2005 = 7990; %million tons of carbon
	%sig0=.13418; % CO2-equivalent emissions-GNP ratio 2005
	p.miu_2005 = 0;
	p.sig0 = p.CO2_2005/p.q0/1000/(1-p.miu_2005);
	p.gsigma=-.158; % Initial growth of sigma per decade
	p.dsig=.00646; % Decline rate of decarbonization per decade
	p.dsig2=0; % Quadratic term in decarbonization
	p.eland0=16; % Carbon emissions from land 2005(GtC per decade)
	p.deland=.2; % Decline rate of land emissions per decade
	
	% Carbon cycle
	p.mat2005=379*2.13; % Concentration in atmosphere 2005 (GtC)
	p.mu2005=1600; % Concentration in upper ocean 2005 (GtC)
	p.ml2005=10010; % Concentration in lower ocean 2005 (GtC)
	p.matPI=278*2.13; % Preindustrial concentration in atmosphere 2005 (GtC)
	p.b=[.88 .12 0 ; .047 .948 .005 ; 0 .001 .999]; % Carbon cycle transition matrix
	
	% Climate model
	p.T2xCO2=3; % Equilibrium temp impact of CO2 doubling oC  
	p.Fex0=.83; % Estimate of 2000 forcings of non-CO2 GHG
	p.Fex1=0.30; % Estimate of 2100 forcings of non-CO2 GHG
	p.forcothmiufactor = 1; % non-CO2 abatement factor (forcoth = forcoth * (1-miu * forcothmiufactor); set to 0 for CO2-only abatement)
	p.forcothlimmiu = 1; % maximum non-CO2 abatement (use <=1 to prevent strong solar radiation management)
	p.Tocean0=.0068; % 2000 lower strat. temp change (C) from 1900
	p.Tatm0=.83; % 2000 atmospheric temp change (C)from 1900 
	p.Trate0=0.03; % 2000 rate of change of temperature (degrees per year)
	p.c1=.208; % Climate-equation coefficient for upper level 
	p.c3=.310; % Transfer coeffic upper to lower stratum
	p.c4=.050; % Transfer coeffic for lower level
	p.FCO22x = 1.66/log2(p.mat2005/p.matPI); % Estimated forcings of equilibrium co2 doubling 
	p.consttcre = 0; % transient climate response to emissions -- should be in units of degree C/Gt C (e.g., 2e-3) if used; note if used non-CO2 forcing is ignored
	
	% Climate damage parameters calibrated for quadratic at 2.5 C for 2105
	p.aa1=0; % damage intercept
	p.aa2=.0028388; % Damage quadratic term
	p.aa3=2.00; % Damage exponent
	p.dam_fcapital=0; % fraction of damages attributed to capital
	p.dam_futility=0; % fraction of damages attributed to effective consumption
	
	% Abatement cost
	p.expcost2=2.8; % Exponent of control cost function 
	p.pback=1.26; % Cost of backstop 2005 000$ per tC 2005   
	p.backrat=2; % Ratio initial to final backstop cost 
	p.gback=.05; % Initial cost decline backstop pc per decade 
	p.limmiu=1; % Upper limit on control rate
	
	% Early retirement additioanl abatement cost
	p.plantlife=1; %lifetime years of energy capital investments (setting to a short value means no early retirement costs)
	p.earlyretcost=.01; % trillion dollars per gigatonne carbon     
	p.abate_futility = 0; % fraction of abatement costs accruing to consumption rather than output   
	
	% Participation
	p.partfract1=1; % Fraction of emissions under control regime 2005 /1      /
	p.partfract2=1; %Fraction of emissions under control regime 2015 /1      /
	p.partfract21=1; % Fraction of emissions under control regime 2205 /1      /
	p.dpartfract=0; % Decline rate of participation                   /0      /
	
	% Availability of fossil fuels
	p.fosslim=6000; %  Maximum cumulative extraction fossil fuels         / 6000  /
	
	% Scaling and inessential parameters
%	p.scale1=194; % Scaling coefficient in the objective function       /194    /
%	p.scale2=381800; % Scaling coefficient in the objective function       /381800 / ;
	%p.scale1=.074;
	%p.scale2=-39370.867;
	p.scale1 = 1; p.scale2 = 0;
	
	% ROI parameters
	p.doROI = 0; % do ROI calculations
	p.ROI_fcons = 1; % fraction of global consumption in ROI
	p.ROI_fpop = 1; % fraction of global population in ROI
	p.ROI_femit = 1; % fraction of global emissions in ROI
	p.ROI_inpartfrac = 1; % is ROI in participation group?
	p.ROI_fmiu = 1; % fraction of abatement borne by ROI
	p.ROI_fliability = 0; % fraction of potential non-ROI libability for which ROI is liable
	p.ROI_damfactor = 1; % damage multiplier for ROI
	p.ROI_limmiu = p.limmiu;
	p.ROI_anchor = [];
	p.ROI_Tstoch = 0; % standard deviation of regional temperature variability
	p.liabilitycap = .9999; % fraction of total transfers possible
	p.liabilityfloor = -Inf; % liability floor; lower to allow reverse liability due to negative cum emissions?
	p.liabilitydamagesref = [];
	
	% environmental goods a la Sterner & Persson 
	p.ecoelasticity = 0.5; % elasticity of substitution
	p.ecoshare = 0; % share of consumption
	
	% new climate module
	p.use_new_climate = 0;
	p.z_mixed = 100; % thickness of mixed layer (m)
	p.z_deep = 3790 - p.z_mixed; % thickness of deep ocean (m)
	p.Cp_mixed = .7 * 3985 * 1030 * p.z_mixed; % heat capacity of surface (J/m^2/K)
	p.sbconst = 5.7e-8 * (pi*1e7); % annualized Stefan-Boltzmann constant (J/m^2/year)
	p.ocean_diffusivity = 1.2e-3 * (pi*1e7); % ocean exchange velocity (m^2/year)
	p.Ndeepboxes = 1;
	p.Tstoch = 0; % standard deviation of regional temperature variability

	% DICE 2007 versions of climate parameters
%	p.Tatm0 = 0.73;
%	p.c1 = 0.220; p.c3 = 0.3;
%	p.mat2005 = 808.0; p.mu2005 = 1255; p.ml2005 = 18365;
%	p.b=[.810712 .189288 0 ; .097213 .852787 .05 ; 0 .003119 .996881];
	
	%%%%%
	
	for i=1:2:length(varargin)
		p.(varargin{i}) = varargin{i+1};
	end
	
	%%%%%

	p.basesavings = 0.2 * ones(size(p.t)); % savings rate to use when not optimizing
	p.lam = p.FCO22x./p.T2xCO2; % climate coefficient

	% population
	%p.Gfacpop = 1 - 1./exp(p.gpop0*((p.t-p.t(1))/10));
	%p.L = p.pop0*(1-p.Gfacpop) + p.Gfacpop*p.popasym;
	p.L = p.pop0 + (p.popasym-p.pop0)*(1 - exp(-p.popadjrate*(p.t-p.t(1))/10));
	
	% productivity
	p.ga = p.ga0 * exp(-p.dela*(p.t-p.t(1)).*exp(-p.deldela*(p.t-p.t(1))));
	p.al = p.a0 * [1 (cumprod((1./(1-p.ga(1:end-1))).^((p.t(2:end)-p.t(1:end-1))/10)))];
	
	% sigma growth
	p.gsig = p.gsigma * exp(-p.dsig*(p.t-p.t(1)) - p.dsig2*(p.t-p.t(1)).^2);
	p.sigma = p.sig0 * [1 (cumprod((1+p.gsig(1:end-1)).^((p.t(2:end)-p.t(1:end-1))/10))) ];
	
	% growth of cost factor
	p.backrat2 = p.backrat/(p.backrat-1);
	p.cost1 = (p.pback*p.sigma/p.expcost2) .* ( (p.backrat2 - 1 + exp(-p.gback .* (p.t-p.t(1))/10))/p.backrat2);
	%p.cost1 = (p.pback*p.sigma/p.expcost2) .* ( (p.backrat - 1 + exp(-p.gback .* (p.t-p.t(1))/10))/p.backrat);
	
	% emissions from deforestation
	p.etree = p.eland0*(1-p.deland).^((p.t-p.t(1))/10);
	
	% average utility social discount factor
	p.rr = 1./((1+p.prstp).^(p.t-p.t(1)));
	
	% exogenous forcing from other greenhouse gases
	p.forcoth = p.Fex0 + (.1*(p.Fex1-p.Fex0)*(p.t-2000)/10).*(p.t<2100) + (p.Fex1-p.Fex0).*(p.t>=2100);

	% participation fraction in the control regime
	p.partfract = ones(size(p.t)) * p.partfract21;
	p.partfract(find(p.t<2200)) = p.partfract21 + (p.partfract2-p.partfract21)*exp(-p.dpartfract*((p.t(find(p.t<2200))-p.t(2))/10));
	p.partfract(1) = p.partfract1;

	p.calcdamages=[];
	
	for i=1:2:length(varargin)
		p.(varargin{i}) = varargin{i+1};
	end
	
	for i=1:2:length(varargin)
		if strcmpi(varargin{i},'calcdamages')
			eval(['p.calcdamages = @(Temperature,Year,Rate,Output,TempHist,Population) ' varargin{i+1} ';']);
		elseif length(varargin{i})>4
			if strcmpi(varargin{i}(1:4),'calc')
				eval(['p.' varargin{i} ' = ' varargin{i+1} ';'])
			end
		end
	end
		
	if length(p.calcdamages)==0
		eval('p.calcdamages = @(Temperature,Year,Rate,Output,TempHist,Population) 1 - 1./(1+p.aa1.*Temperature+p.aa2.*abs(Temperature).^p.aa3);');
	end

	if length(p.ROI_fpop)==1
		p.ROI_fpop = repmat(p.ROI_fpop,1,length(p.t));
	end
	if length(p.ROI_fcons)==1
		p.ROI_fcons = repmat(p.ROI_fcons,1,length(p.t));
	end
	if length(p.ROI_femit)==1
		p.ROI_femit = repmat(p.ROI_femit,1,length(p.t));
	end
	if length(p.ROI_fliability)==1
		p.ROI_fliability = repmat(p.ROI_fliability,1,length(p.t));
	end

	p.ROI_partfrac = min(1,p.partfract./p.ROI_femit);

	if length(p.ROI_fmiu)==1
		p.ROI_fmiu = ((p.ROI_femit.*p.ROI_inpartfrac)>=p.partfract) + ((p.ROI_femit.*p.ROI_inpartfrac)<p.partfract) .* ((p.ROI_femit.*p.ROI_inpartfrac)./p.partfract);
	end

	if p.ecoelasticity == 1
		p.ecoelasticity = 0.9999;
	end
	if length(p.ecoshare)>1
		p.ecoshare = repmat(p.ecoshare,1,length(p.t));
	end
	if length(p.ecoelasticity)>1
		p.ecoelasticity = repmat(p.ecoelasticity,1,length(t));
	end

	if length(p.liabilitydamagesref)==0
		if length(p.T2xCO2)>1
			for i=1:length(p.t)
			 	p.liabilitydamagesref(:,i) = p.calcdamages(3 - (3-p.Tatm0) * (exp(-(p.t(i)-p.t(1))/30)));
			 end
		else
			 	p.liabilitydamagesref = p.calcdamages(3 - (3-p.Tatm0) .* exp(-(p.t-p.t(1))/30 ));
		end
	 end
	 
	 if length(p.ROI_anchor) == 0
	 	p.ROI_anchor = zeros(2,length(p.t));
	end
	
	if p.ROI_Tstoch>0
		p.ROI_Tstoch_seeds = p.ROI_Tstoch .* randn(length(p.T2xCO2),length(p.t));
	end
	
	if p.Tstoch>0
		p.Tstoch_seeds = p.Tstoch .* randn(length(p.T2xCO2),length(p.t));
	end