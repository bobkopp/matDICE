classdef DICEParameters < handle
% DICE Parameters class
%
% Last updated by Bob Kopp rkopp-at-alumni.caltech.edu, 7 April 2012
	
	properties (SetAccess = private)
		t = 2005:10:2305; % years

		% Discounting
		elasmu = 1.5; %Elasticity of marginal utility of consumption
		prstp=0.015; %Initial rate of social time preference per year
	
		basesavings = 0.2;
	
		% Population and technology
		pop0=6411; % 2005 world population millions
		popasym=8700; % Asymptotic population
		popadjrate=0.59; % Population adjustment rate per decade
		al; % Total factor productivity
		ga0=.160; % Initial growth rate of productivity per decade
		dela=.0094; % Decline rate of productivity growth per decade
		deldela=0.0019; % Rate of decline of decline rate in productivity growth
		dk=.100; % Depreciation rate on capital per year
		gama=.300; % Capital elasticity in production function
		q0=57.09; % 2005 world gross output trill 2005 US international dollars % matches world bank figures
		k0=97.3/55.3 * 57.09; %2005 value capital trill 2005 US dollars
		minoutpc = 730; % dollars per person per year
		mincons = 730; % dollars per person per year
		popfeedback = 0; % population feedback on or off (do we react to output less than the minimum by raising output or cutting population?)
		allowrecoveryfromcrash = 0; % are crashes to minimum output forever?
		L; % population
		
		% Emissions
		
		CO2_2005 = 7990; %million tons of carbon
		miu_2005 = 0;
		gsigma=-.158; % Initial growth of sigma per decade
		dsig=.00646; % Decline rate of decarbonization per decade
		dsig2=0; % Quadratic term in decarbonization
		eland0=16; % Carbon emissions from land 2005(GtC per decade)
		deland=.2; % Decline rate of land emissions per decade
		sigma; % carbon intensity
		etree; %emissions from deforestation
		
		% Carbon cycle
		mat2005=379*2.13; % Concentration in atmosphere 2005 (GtC)
		mu2005=1600; % Concentration in upper ocean 2005 (GtC)
		ml2005=10010; % Concentration in lower ocean 2005 (GtC)
		matPI=278*2.13; % Preindustrial concentration in atmosphere 2005 (GtC)
		b=[.88 .12 0 ; .047 .948 .005 ; 0 .001 .999]; % Carbon cycle transition matrix
		
		% Climate model
		T2xCO2=3; % Equilibrium temp impact of CO2 doubling oC  
		FCO22x=1.66/log2(379/278); % Estimated forcings of equilibrium co2 doubling 	
		Fex0=.83; % Estimate of 2000 forcings of non-CO2 GHG
		Fex1=0.30; % Estimate of 2100 forcings of non-CO2 GHG
		forcothmiufactor = 1; % non-CO2 abatement factor (forcoth = forcoth * (1-miu * forcothmiufactor); set to 0 for CO2-only abatement)
		forcothlimmiu = 1; % maximum non-CO2 abatement (use <=1 to prevent strong solar radiation management)
		forcoth; % non-CO2 forcing
		Tocean0=.0068; % 2000 lower strat. temp change (C) from 1900
		Tatm0=.83; % 2000 atmospheric temp change (C)from 1900 
		Trate0=0.03; % 2000 rate of change of temperature (degrees per year)
		c1=.208; % Climate-equation coefficient for upper level 
		c3=.310; % Transfer coeffic upper to lower stratum
		c4=.050; % Transfer coeffic for lower level
		
		% Climate damage parameters calibrated for quadratic at 2.5 C for 2105
		aa1=0; % damage intercept
		aa2=.0028388; % Damage quadratic term
		aa3=2.00; % Damage exponent
		dam_fcapital=0; % fraction of damages attributed to capital
		dam_futility=0; % fraction of damages attributed to effective consumption
		
		% Abatement cost
		expcost2=2.8; % Exponent of control cost function 
		pback=1.26; % Cost of backstop 2005 000$ per tC 2005   
		backrat=2; % Ratio initial to final backstop cost 
		gback=.05; % Initial cost decline backstop pc per decade 
		limmiu=1; % Upper limit on control rate
		
		% Early retirement additioanl abatement cost
		plantlife=1; %lifetime years of energy capital investments (setting to a short value means no early retirement costs)
		earlyretcost=.01; % trillion dollars per tonne carbon        
		
		% Participation
		partfract1=1; % Fraction of emissions under control regime 2005 /1      /
		partfract2=1; %Fraction of emissions under control regime 2015 /1      /
		partfract21=1; % Fraction of emissions under control regime 2205 /1      /
		dpartfract=0; % Decline rate of participation                   /0      /
		
		% Availability of fossil fuels
		fosslim=6000; %  Maximum cumulative extraction fossil fuels         / 6000  /
		
		% Scaling and inessential parameters
	%	scale1=194; % Scaling coefficient in the objective function       /194    /
	%	scale2=381800; % Scaling coefficient in the objective function       /381800 / ;
		%scale1=.074;
		%scale2=-39370.867;
		scale1 = 1; scale2 = 0;
		ROI_fcons = 1; % fraction of global consumption in ROI
		ROI_fpop = 1; % fraction of global population in ROI 
		
		% environmental goods a la Sterner & Persson 
		ecoelasticity = 0.5; % elasticity of substitution
		ecoshare = 0; % share of consumption

	
		calcdamagesstr = '@(Temperature,Year,Rate,Output,TempHist,Population) 1 - 1./(1+p.aa1.*Temperature+p.aa2.*abs(Temperature).^p.aa3);';
		
		spec = {};
	end
	
	properties (Dependent = true, SetAccess = private)
		lam; % climate coefficient
%		Gfacpop; % population growth factor
		sig0; % CO2-equivalent emissions-GNP ratio 2005
		a0; % Initial level of total factor productivity
		ga; % productivity growth factor
		gsig; % sigma growth factor
		cost1; % abatement cost factor
		rr; % time discount factor
		partfract; % participation fraction
	end
	
	methods
	
		function obj=DICEParameters(varargin)
			updateL(obj);
			updateal(obj);
			updatesigma(obj);
			updateetree(obj);
			updateforcoth(obj);		

			if nargin>0
				obj.setPropValue(varargin{:});
			end

			updateL(obj);
			updateal(obj);
			updatesigma(obj);
			updateetree(obj);
			updateforcoth(obj);		

			if nargin>0
				obj.setPropValue(varargin{:});
			end

		end
		
		function o=setPropValue(obj,varargin)
			for i=2:2:length(varargin)
				if isprop(obj,varargin{i-1})
					obj.(varargin{i-1})=varargin{i};
					obj.spec={obj.spec{:},varargin{i-1},varargin{i}};
				else
					disp(['No property ' varargin{i-1}]);
				end
			end

		end
		
		function u=get.partfract(obj)
			u = ones(size(obj.t)) * obj.partfract21;
			u(find(obj.t<2200)) = obj.partfract21 + (obj.partfract2-obj.partfract21)*exp(-obj.dpartfract*((obj.t(find(obj.t<2200))-obj.t(2))/10));
			u(1) = obj.partfract1;
		end

		function u=updateL(obj)
			obj.L = obj.pop0 + (obj.popasym-obj.pop0)*(1 - exp(-obj.popadjrate*(obj.t-obj.t(1))/10));
		end
		
		function u=updateal(obj)
			obj.al = obj.a0 * [1 cumprod(1./(1-obj.ga(1:end-1)))];
		end

		function u = updatesigma(obj)
				obj.sigma = obj.sig0 * [1 cumprod(1+obj.gsig(1:end-1)) ];
		end

		function u = updateetree(obj)
			obj.etree = obj.eland0*(1-obj.deland).^((obj.t-obj.t(1))/10);
		end
				
		function u = updateforcoth(obj)
			obj.forcoth = obj.Fex0 + (.1*(obj.Fex1-obj.Fex0)*(obj.t-2000)/10).*(obj.t<2100) + (obj.Fex1-obj.Fex0).*(obj.t>=2100);
		end
		
		function u = set.a0(obj,value)
		end
		
		function u = set.sig0(obj,value)
		end
		
		function u = set.sigma(obj,value)
			obj.sigma = value;
			obj.CO2_2005 = obj.sigma(1) * obj.q0 * 1000 * (1-obj.miu_2005);
		end
		
		function u = get.a0(obj)
			u = obj.q0/(obj.k0^obj.gama * obj.pop0^(1-obj.gama));
		end

		function u = get.sig0(obj)
			u = obj.CO2_2005/obj.q0/1000/(1-obj.miu_2005);
		end
		
		function calclam = get.lam(obj)
			calclam = obj.FCO22x./obj.T2xCO2;
		end
		
		function u = get.basesavings(obj)
			if length(obj.basesavings)==1
				u = obj.basesavings * ones(size(obj.t));
			else
				u = obj.basesavings;
			end
		end

		function u = get.ROI_fpop(obj)
			if length(obj.ROI_fpop)==1
				u = obj.ROI_fpop * ones(size(obj.t));
			else
				u = obj.ROI_fpop;
			end
		end
	
		function u = get.ROI_fcons(obj)
			if length(obj.ROI_fcons)==1
				u = obj.ROI_fcons * ones(size(obj.t));
			else
				u = obj.ROI_fcons;
			end
		end
	
		function u = get.ecoshare(obj)
			if length(obj.ecoshare)>1
				u = repmat(obj.ecoshare,1,length(obj.t));
			else
				u = obj.ecoshare;
			end
		end

		function u = get.ecoelasticity(obj)
			u = obj.ecoelasticity - 1e-4*(obj.ecoelasticity==1);
			if length(u)>1
				u = repmat(u,1,length(obj.t));
			end
		end
		
%		function calcGfacpop = get.Gfacpop(obj)
%			calcGfacpop = 1 - 1./exp(obj.gpop0*((obj.t-obj.t(1))/10));
%		end
		
		function calcga = get.ga(obj)
			calcga = obj.ga0 * exp(-obj.dela*(obj.t-obj.t(1)).*exp(-obj.deldela*(obj.t-obj.t(1))));
		end
		
		function calcgsig = get.gsig(obj)
			calcgsig = obj.gsigma * exp(-obj.dsig*(obj.t-obj.t(1)) - obj.dsig2*(obj.t-obj.t(1)).^2);
		end
		
		function calccost1 = get.cost1(obj)
			calccost1 = (obj.pback*obj.sigma/obj.expcost2) .* ( (obj.backrat - 1 + exp(-obj.gback .* (obj.t-obj.t(1))/10))/obj.backrat);
		end
		
		function calcrr = get.rr(obj)
			calcrr = 1./((1+obj.prstp).^(obj.t-obj.t(1)));
		end		
		
		function u=calcdamages(p,varargin)
			eval(['damfunc = ' p.calcdamagesstr]);
			u=damfunc(varargin{:});
		end		
	end
end