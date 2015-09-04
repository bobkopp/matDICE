function p=DICEParameters(varargin)

% p=DICEParameters([parameter],[value],...)
%
% Example:
%     p = DICEParameters('elasmu',0,'prstp',0.03);
%
% Default values largely match DICE-2013R.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Sep 03 21:27:42 EDT 2015
%
% Copyright (C) 2015 by Robert E. Kopp; distributed under GNU GPL v3

    p.spec = varargin;

    p.t=2010:5:2400;
    
    % Preferences
    
    p.elasmu=1.45; %Elasticity of marginal utility of consumption
    p.prstp=0.015; %Initial rate of social time preference per year
    
    % Population and technology
    p.gama=.300; % Capital elasticity in production function
    p.pop0=6838; % Initial world population (millions)
    p.popadj = 0.134; % Growth rate to calibrate to 2050 pop projection
    p.popasym=10500; % Asymptotic population (millions)
                     %p.gpop0=.115; % Growth rate of population per decade
                     %p.popadjrate=0.59; % Population adjustment rate per decade
    p.dk=.100; % Depreciation rate on capital per year
    p.q0=63.69; % Initial world gross output (trill 2005 USD)
    p.k0=135; % Initial capital value (trill 2005 USD)           
    p.a0=3.80; % Initial level of total factor productivity
               %p.a0 = p.q0/(p.k0^p.gama * p.pop0^(1-p.gama));
    p.ga0=.079; % Initial growth rate for TFP per 5 years
    p.dela=.006; % Decline rate of TFP Per 5 years
    
    %p.deldela=0.0019; % Rate of decline of decline rate in productivity growth
    p.minoutpc = 730; % dollars per person per year
    p.mincons = 730; % dollars per person per year
    p.popfeedback = 0; % population feedback on or off (do we react to output less than the minimum by raising output or cutting population?)
    p.allowrecoveryfromcrash = 0; % are crashes to minimum output forever?

    % Emissions
    
    p.gsigma=-.01; % Initial growth of sigma (per year)
    p.dsig=-0.001; % Decline rate of decarbonization (per period)
    p.eland0=3.3; % Carbon emissions from land 2010 (GtCO2 per year)
    p.deland=.2; % Decline rate of land emissions (per period)
    p.e0 = 33.61; % Industrial emissions 2010 (GtCO2 per year);
    p.miu0 = 0.039; % Initial emissions control rate for base case 2010
    
    % Carbon cycle
    p.mat0=389*2.13; % Initial Concentration in atmosphere 2010 (GtC)
    p.mu0=1527; % Concentration in upper ocean 2010 (GtC)
    p.ml0=10010; % Concentration in deep ocean 2010 (GtC)
    p.mateq=278*2.13; % Equilibrium concentration atmosphere  (GtC)
    p.mueq = 1350; % Equilibrium concentration in upper ocean (GtC)
    p.mleq = 10000; % Equilibrium concentration in deep ocean (GtC)
    
    % Carbon cycle transition matrix
    p.b12=0.088; p.b23 = 0.00250;
    b21 = p.b12*p.mateq/p.mueq;
    b32 = p.b23*p.mueq/p.mleq;
       
    % Climate model
    p.T2xCO2=2.9; % Equilibrium temp impact (oC per doubling CO2)  
    p.Fex0=0.25; % 2010 forcings of non-CO2 GHG (Wm-2)
    p.Fex1=0.70; % 2100 forcings of non-CO2 GHG (Wm-2)
    
    p.forcothmiufactor = 1; % non-CO2 abatement factor (forcoth = forcoth * (1-miu * forcothmiufactor); set to 0 for CO2-only abatement)
    p.forcothlimmiu = 1; % maximum non-CO2 abatement (use <=1 to prevent strong solar radiation management)
    
    p.Tocean0=.0068; % Initial lower stratum temp change (C from 1900)
    p.Tatm0=0.80; % Initial atmospheric temp change (C from 1900) 
    p.Trate0=0.03; % 2000 rate of change of temperature (degrees per year)

    p.c10 = 0.098; %Initial climate equation coefficient for upper level
    p.c1beta = 0.01243; % Regression slope coefficient(SoA~Equil TSC)
                        % NOTE THAT THIS APPEARS TO BE OF THE WRONG SIGN IN DICE-2013R
                        % C1 IS THE INVERSE OF OCEAN HEAT CAPACITY, SO A POSITIVE C1BETA
                        % MEANS THAT INCREASING CLIMATE SENSITIVITY DECREASES OCEAN HEAT CAPACITY
    
    p.c3=.088; % Transfer coefficient upper to deep ocean
    p.c4=.025; % Transfer coefficient for lower level
    p.c4toc3=p.c4/p.c3; % ratio of deep to shallow heat capacity (this seems quite large in DICE default)
    
    p.FCO22x = 3.8; % Forcings of equilibrium CO2 doubling (Wm-2)       
    
    p.consttcre = 0; % transient climate response to emissions -- should be in units of degree C/Gt C (e.g., 2e-3) if used; note if used non-CO2 forcing is ignored
    
    % Climate damage parameters calibrated for quadratic at 2.5 C for 2105
    p.aa1=0; % damage intercept
    p.aa2=.00267; % Damage quadratic term
    p.aa3=2.00; % Damage exponent
    p.dam_fcapital=0; % fraction of damages attributed to capital
    p.dam_futility=0; % fraction of damages attributed to effective consumption
    
    % Abatement cost
    p.expcost2=2.8; % Exponent of control cost function 
    p.pback=344; % Cost of backstop 2005 $ per tCO2 2010   
    p.gback=.025; % Initial cost decline backstop pc per period 
    p.limmiu=1.2; % Upper limit on control rate after 2150
    p.backrat=2; % Ratio initial to final backstop cost 
    p.tnopol=2230; %    Period before which no emissions controls base  / 45   /
    p.cprice0=1.0; %   Initial base carbon price (2005$ per tCO2)      / 1.0  /
    p.gcprice=0.02; %   Growth rate of base carbon price per year       /.02   /


    
    % Early retirement additional abatement cost
    p.plantlife=1; %lifetime years of energy capital investments (setting to a short value means no early retirement costs)
    p.earlyretcost=.01; % trillion dollars per gigatonne carbon     
    p.abate_futility = 0; % fraction of abatement costs accruing to consumption rather than output   
    
    % Participation
    p.periodfullpart=2110; % Period at which have full participation           /21  /
    p.partfract2010=1; %  Fraction of emissions under control in 2010       / 1  /
    p.partfractfull=1; % Fraction of emissions under control at full time  / 1  /
    
    % Availability of fossil fuels
    p.fosslim=6000; %  Maximum cumulative extraction fossil fuels (GtC)        / 6000  /
    
    % Scaling and inessential parameters
    %	p.scale1=194; % Scaling coefficient in the objective function       /194    /
    %	p.scale2=381800; % Scaling coefficient in the objective function       /381800 / ;
    %p.scale1=0.016408662;
    %p.scale2=-3855.106895;
    p.scale1 = 1; p.scale2 = 0;
    
    % Region Of Interest (ROI) parameters (Kopp & Mignone 2013)
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

    % for climate experiments
    p.overrideForcing=[];
    p.overrideCO2=[];
    p.overrideCarbon=[];
    
    % DICE 2007 versions of climate parameters
    %	p.Tatm0 = 0.73;
    %	p.c1 = 0.220; p.c3 = 0.3;
    %	p.mat2005 = 808.0; p.mu2005 = 1255; p.ml2005 = 18365;
    %	p.b=[.810712 .189288 0 ; .097213 .852787 .05 ; 0 .003119 .996881];
    
    %%%%%
    
    p.optlrsav = (p.dk + .004)/(p.dk + .004*p.elasmu + p.prstp)*p.gama;

    for i=1:2:length(varargin)
        p.(varargin{i}) = varargin{i+1};
    end

    %%%%%

    
    p.b=[(1-p.b12) p.b12 0 ; b21 1-b21-p.b23 p.b23 ; 0 b32 1-b32]; 
    p.sig0 = p.e0/(p.q0-(1-p.miu0));
    p.c4=p.c4toc3.*p.c3; % preserve ratio of deep to shallow heat capacity


    %%%%
    
    %l(t)          Level of population and labor
        %al(t)         Level of total factor productivity
        %sigma(t)      CO2-equivalent-emissions output ratio
        %rr(t)         Average utility social discount rate
        %ga(t)         Growth rate of productivity from
        %forcoth(t)    Exogenous forcing for other greenhouse gases
            %gl(t)         Growth rate of labor
        %gcost1        Growth of cost factor
        %gsig(t)       Change in sigma (cumulative improvement of energy efficiency)
        %etree(t)      Emissions from deforestation
        %cost1(t)      Adjusted cost for backstop
            %partfract(t)  Fraction of emissions in control regime
        %lam           Climate model parameter
        %gfacpop(t)    Growth factor population
        %pbacktime(t)  Backstop price
        %optlrsav      Optimal long-run savings rate used for transversality
        %cpricebase(t) Carbon price in base case
        
        %%%%

        p.lam = p.FCO22x./p.T2xCO2; % climate feedback coefficient
        p.c1 =  p.c10 + p.c1beta*(p.T2xCO2-2.9);
 
    % population - change relative to DICE for continuous time
    q=p.popasym./p.pop0; alf = p.popadj;
    stepsize=p.t(2)-p.t(1);
    u=-(1/stepsize)*log((q.^alf-q)/(1-q));
    p.poptau=1./u;
    popdelta=(p.popasym-p.pop0).*exp(-(p.t-p.t(1))/p.poptau);
    p.L=p.popasym-popdelta;
    
    % productivity
    p.ga = p.ga0 * exp(-p.dela*(p.t-p.t(1)));
    p.al = p.a0 * [1 cumprod((1./(1-p.ga(1:end-1))))];
    
    % sigma growth
    p.gsig = p.gsigma * (1+p.dsig).^(p.t-p.t(1));
    p.sigma = p.sig0 * [1 cumprod(exp(p.gsig(1:end-1)*stepsize)) ];
    
    % growth of cost factor
    p.pbacktime = p.pback*(1-p.gback).^((p.t-p.t(1))/stepsize);
    p.cost1 = p.pbacktime.*p.sigma/p.expcost2/1000;
 
 
    % emissions from deforestation
    p.etree = p.eland0*(1-p.deland).^((p.t-p.t(1))/stepsize);
    
    % average utility social discount factor
    p.rr = 1./((1+p.prstp).^(p.t-p.t(1)));
    
    % exogenous forcing from other greenhouse gases
    p.forcoth = p.Fex0 + ((1/18)*(p.Fex1-p.Fex0)*(p.t-p.t(1))/stepsize).*(p.t<2100) + (p.Fex1-p.Fex0).*(p.t>=2100);

    p.basesavings = 0.2 * ones(size(p.t)); % savings rate to use when not optimizing


    % participation fraction in the control regime
    p.partfract = ones(size(p.t)) * p.partfractfull;
    p.partfract(p.t<p.periodfullpart) = p.partfract2010 + (p.partfractfull-p.partfract2010)*(p.t(p.t<p.periodfullpart)-2010)./(p.periodfullpart-2010);

    p.cpricebase = p.cprice0+(1+p.gcprice).^(p.t-p.t(1));
    
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