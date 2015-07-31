function p=RICEParameters(varargin)

% STILL NEED TO DO NEGISHI WEIGHTS

	p.t=2005:10:2305;

	% US EU Japan Russia Eurasia China India MidEast Africa LA OHI OthAsia
	
	% Preferences
	
	p.elasmu=1.5; %Elasticity of marginal utility of consumption  /2.0/   
	p.prstp=0.015; %Initial rate of social time preference per year /0.015/
	
	% Population and technology
	p.pop0=[297 490 128 143 156 1305 1095 413 764 555 129 937];
	p.gpop0=[0.093 0.043 -0.013 -0.037 0.011 0.062 0.135 0.190 0.215 0.106 0.073 0.143];
	p.popasym=[478 521 93 105 137 1298 1399 835 1573 685 150 1425];

	% Productivity/Growth

	%p.ga0=.230; % Initial growth rate for technology per decade   /.092      /
	p.growth0 = [0.0151 0.0163 0.0138 0.0260 0.0247 0.0714 0.0455 0.0235 0.0365 0.0293 0.0188 0.0279]% Initial growth rate;
	p.dela=.1; % Decline rate of TFP per decade       /.1     /
	p.tfpconv = 0.10; % TFP convergence rate (per decade) /.1/
	p.glongrun = 0.0033; %long run growth rate (per year) /.0.0033/
	p.convratio = [1 .9 .9 .6 .6 .6 .5 .5 .4 .7 .9 .6]; %GDP/cap convergence ratio
	
	% Output and production
	
	p.dk=.100; % Depreciation rate on capital per year           /.100     /
	p.gama=.300; % Capital elasticity in production function       /.300     /
	p.q0=[ 12.398 	 13.031 	 3.870 	 1.698 	 0.807 	 5.333 	 2.441 	 3.480 	 1.301 	 4.558 	 3.842 	 2.619 ] % 2005 world gross output trill 2005 US dollars
	p.k0=[ 22.851 	 23.302 	 7.133 	 2.787 	 1.362 	 9.261 	 4.119 	 5.416 	 2.135 	 7.693 	 6.870 	 4.420 ]; %2005 value capital trill 2005 US dollars        /137.     /
	p.a0 = p.q0./(p.k0.^p.gama .* p.pop0.^(1-p.gama)); % Initial level of total factor productivity
	
	% Emissions
	
	p.CO2_2005 = [ 1662 	 1146 	 370 	 431 	 257 	 1601 	 405 	 590 	 191 	 412 	 542 	 364 ];
	%p.miu_2005 = repmat(0.005,1,length(p.CO2_2005));
	%p.sig0 = p.CO2_2005./p.q0/1000./(1-p.miu_2005);
	p.sig0 = [0.134065723 0.08794763 0.095591421 0.253973604 0.317917935 0.300144622 0.165989167 0.169441796 0.147153133 0.090440835 0.141012818 0.138956043];
	p.gsig0 = [-0.0175 -0.0145 -0.0175 -0.0173 -0.0240 -0.0249 -0.0221 -0.0176 -0.0221 -0.0150 -0.0187 -0.0190];
	p.gsigaddfactor = [0.903 0.771 1.076 0.882 0.854 1.017 1.181 0.985 0.842 0.962 0.956 1.314]; %2015 add factors
	p.gsigma=-.0025; % Decline rate sigma growth              /-.0025    /
	p.dsig=.1; % Decline rate of decarbonization per decade      /.1   /
	p.eland0=[0 0 0 0 0 0 0 0 0.3 0.6 0 0.7] * 10; % Carbon emissions from land 2005(GtC per decade) / 11.000  /
	p.deland=.2; % Decline rate of land emissions per decade
	
	% Carbon cycle
	p.mat2005=379*2.13; % Concentration in atmosphere 2005 (GtC)          /808.9   /
	p.mu2005=1600; % Concentration in upper ocean 2005 (GtC)        /1255     /
	p.ml2005=10010; % Concentration in lower ocean 2005 (GtC)        /18365    /
	p.matPI=278*2.13; % Preindustrial concentration in atmosphere 2005 (GtC) /596.4/
	%p.b=[.810712 .189288 0 ; .097213 .852787 .05 ; 0 .003119 .996881]; % Carbon cycle transition matrix
	p.b=[.88 .12 0 ; .047 .948 .005 ; 0 .001 .999];
	
	% Climate model
	p.T2xCO2=3.2; % Equilibrium temp impact of CO2 doubling oC      / 3.2 /
	p.Fex0=-.06; % Estimate of 2000 forcings of non-CO2 GHG        / -.06   /
	p.Fex1=0.30; % Estimate of 2100 forcings of non-CO2 GHG        / 0.30   /
	p.Tocean0=.0068; % 2000 lower strat. temp change (C) from 1900     /.0068   /
	p.Tatm0=.83; % 2005 atmospheric temp change (C)from 1900       /.83   /
	p.c1=.190; % Climate-equation coefficient for upper level    /.220    /
	p.c3=.310; % Transfer coeffic upper to lower stratum         /.300    /
	p.c4=.050; % Transfer coeffic for lower level                /.050    /
	%p.FCO22x=3.8; % Estimated forcings of equilibrium co2 doubling  /3.8     /
	p.FCO22x = 1.66/log2(p.mat2005/p.matPI);
	
	% Climate damage parameters calibrated for quadratic at 2.5 C for 2105
	p.aa1=0; % damage intercept                                / 0.00000    /
	p.aa2=.0028388; % Damage quadratic term                           /  0.0028388 /
	p.aa3=2.00; % Damage exponent                                 / 2.00       /
	
	% Abatement cost
	p.expcost2=2.8; % Exponent of control cost function               /2.8   /
	p.pback=1.26; % Cost of backstop 2005 000$ per tC 2005          /1.17  /
	p.backyear = 2250; %backstop inflection year /2250/
	p.backrat=0.1; % Ratio asymptotic to initial backstop price            / 0.1    /
	p.gback=.05; % Initial cost decline backstop pc per decade     /.05   /
	p.limmiu=1; % Upper limit on control rate                     / 1    /
	p.backstopreg = [0.9 1.4 1.4 0.6 0.6 0.7 1.1 1.0 1.1 1.3 1.1 1.2]; % ratio of regional backstop price to world
	
	% Participation
%	p.partfract1=1; % Fraction of emissions under control regime 2005 /1      /
%	p.partfract2=1; %Fraction of emissions under control regime 2015 /1      /
%	p.partfract21=1; % Fraction of emissions under control regime 2205 /1      /
%	p.dpartfract=0; % Decline rate of participation                   /0      /
	
	% Availability of fossil fuels
	p.fosslim=6000; %  Maximum cumulative extraction fossil fuels         / 6000  /
	
	% Scaling and inessential parameters
%	p.scale1=194; % Scaling coefficient in the objective function       /194    /
%	p.scale2=381800; % Scaling coefficient in the objective function       /381800 / ;
	%p.scale1=.074;
	%p.scale2=-39370.867;
	p.scale1 = 1; p.scale2 = 0;
	
	% environmental goods a la Sterner & Persson 
	p.ecoelasticity = .5; % elasticity of substitution
	p.ecoshare = 0; % share of consumption
	p.ecodamagec1 = 0; % scalar in eco damage function
	p.ecodamagec2 = 2; % power in eco damage function

	%%%%%
	
	for i=1:2:length(varargin)
		p.(varargin{i}) = varargin{i+1};
	end
	
	%%%%%

	p.lam = p.FCO22x/p.T2xCO2; % climate coefficient

	% population
	p.Gfacpop = 1 - 1./exp(p.gpop0'*((p.t-p.t(1))/10));
	p.L = bsxfun(@times,p.pop0',(1-p.Gfacpop)) + bsxfun(@times,p.popasym',p.Gfacpop);
	
	% sigma growth
	p.gsig = p.gsigma + (p.gsig0-p.gsigma)' * (1-p.dsig).^((p.t(2:end)-p.t(2))/10);
	p.sigma = bsxfun(@times,p.sig0',[ones(size(p.q0')) cumprod(exp(bsxfun(@times,p.gsig,p.t(2:end)-p.t(1:end-1))),2)]);
	p.sigma(:,2:end) = bsxfun(@times,p.gsigaddfactor',p.sigma(:,2:end));
	
	% productivity

	p.q0percap = p.q0./p.pop0*1000;
	p.Gqpercap(1,:) = [exp((p.t(2:end)-p.t(1:end-1))*p.glongrun).*exp((p.t(2:end)-p.t(1:end-1)).*(p.growth0(1)-p.glongrun).*exp(-p.dela*(p.t(2:end)-p.t(2))/10))];
	p.qpercap(1,:) = p.q0percap(1) * [1 cumprod(p.Gqpercap(1,:))];
	
	p.Gqpercap(2:length(p.q0),1) = [exp((p.t(2)-p.t(1))*p.growth0(2:end)')];
	p.qpercap(2:length(p.q0),1:2) = bsxfun(@times,p.q0percap(2:end)',[ones(length(p.q0)-1,1) p.Gqpercap(2:end,1)]);
	for i=3:length(p.t)
		p.qpercap(2:length(p.q0),i) = p.qpercap(2:end,i-1)*p.Gqpercap(1,i-1).*(p.qpercap(1,i-1)./p.qpercap(2:end,i-1).*p.convratio(2:end)').^(p.tfpconv*(p.t(i)-p.t(i-1))/10);
	end
	p.Gqpercap=p.qpercap(:,2:end)./p.qpercap(:,1:end-1);
	p.tfpgrowth = log(p.Gqpercap).*(1-p.gama);
	p.al = bsxfun(@times,p.a0',cumprod([ones(length(p.q0),1) exp(p.tfpgrowth)],2));
	
	p.q = p.qpercap/1000.*p.L;
	
	% growth of cost factor
	p.backcost= (p.t<p.backyear) .* p.pback.*(p.backrat+(1-p.backrat).*(1-p.gback).^(.1*(p.t-p.t(1))));
	p.backcost=p.backstopreg'*p.backcost;
	p.cost1 = (p.backcost.*p.sigma/p.expcost2);
	% emissions from deforestation
	p.etree = bsxfun(@times,p.eland0',(1-p.deland).^((p.t-p.t(1))/10));
	
	% average utility social discount factor
	p.rr = 1./((1+p.prstp).^(p.t-p.t(1)));
	
	% exogenous forcing from other greenhouse gases
	p.forcoth = p.Fex0 + (.1*(p.Fex1-p.Fex0)*(p.t-2000)/10).*(p.t<2100) + (p.Fex1-p.Fex0).*(p.t>=2100);
	
	
end