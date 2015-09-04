function [Base,OptimAbatement,FullOptim] = matDICE(varargin)

% [Base,OptimAbatement,FullOptim] = matDICE(['parameter',value,'parameter',value,...])
%
% FullOptim is calculated allowing optimization of both abatement and savings rate.
% OptimAbatement is calculated allowing optimization only of abatement.
% Base is caulculated with both fixed.
%
% The switch 'doMonotonicity' requires abatement to monotonically increase.
% The parameter 'FixedAbatement' locks abatement for the first n periods at the values
% specified by the following vector.
% Other parameters are as in DICEParameters.
%
% Example: Calculation at 3% fixed discount rate
%
%         [Base,OptimAbatement,FullOptim] = matDICE('elasmu',0,'prstp',.03);
%         [BaseNoDamages] = matDICE('elasmu',0,'prstp',.03,'aa2',0);
%
%         clf;
%         subplot(2,1,1); hold on;
%         plot(Base.times,Base.ppmCO2,'r');
%         plot(FullOptim.times,FullOptim.ppmCO2,'b');
%         plot(OptimAbatement.times,OptimAbatement.ppmCO2,'g');
%         title('CO_2'); ylabel('ppm CO_2');
%
%         subplot(2,1,2); hold on;
%         plot(Base.times,Base.Output,'r');
%         plot(FullOptim.times,FullOptim.Output,'b');
%         plot(OptimAbatement.times,OptimAbatement.Output,'g');
%         plot(BaseNoDamages.times,BaseNoDamages.Output,'k');
%         title('GDP'); ylabel('Output (US$ trillion)');
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Sep 03 21:46:21 EDT 2015
% Copyright (C) 2015 by Robert E. Kopp; distributed under GNU GPL v3
	
	doMonotonicity=0;
	setFixedAbatement = 0;
	coarseners=0;
	baseabate = [];
	basesave = [];
	MaxFunEvals = 8000;
	incrementalemissions=[];
	startabate =  [];
	dispmode = 'final';
	useparallel = 'always';
	usenlopt = 0;
	TolX = 1e-4;
	TolFun = 1e-6;
	algorithm = 'sqp';
	
	i=1; j=1;
	spec={};
	while(i<=length(varargin))
		if strcmpi(varargin{i},'monotonicAbatement');
			doMonotonicity=1;
			i=i+1;
		elseif strcmpi(varargin{i},'FixedAbatement');
			FixedAbatement=varargin{i+1};
			setFixedAbatement = 1;
			i=i+2;
		elseif strcmpi(varargin{i},'BaseAbatement');
			baseabate = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'BaseSavings');
			basesave = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'incrementalemissions');
			incrementalemissions = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'startabatement');
			startabate = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'Display');
			dispmode = varargin{i+1};
                        i=i+2;
                elseif strcmpi(varargin{i},'coarseners');
                    coarseners = varargin{i+1};
                    coarseners = sort(coarseners,1,'descend');
			i=i+2;
		elseif strcmpi(varargin{i},'MaxFunEvals');
			MaxFunEvals = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'TolX');
			TolX = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'TolFun');
			TolFun = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'algorithm');
			algorithm = varargin{i+1};
			i=i+2;
		elseif strcmpi(varargin{i},'nlopt')
			usenlopt = 1;
			i=i+1;
		else
			spec{j} = varargin{i};
			spec{j+1} = varargin{i+1};
			%eval(['p.' varargin{i} '=' num2str(varargin{i+1}) ';']);
			i=i+2;
			j=j+2;
		end
	end
	p=DICEParameters(spec{:});
	timesteps=length(p.t);
	if ~setFixedAbatement
		FixedAbatement = p.miu0;
	end

	if length(basesave)==0
		basesavings = p.basesavings;
	elseif length(basesave)==1
		basesavings = basesave*ones(size(p.t));
	else
		basesavings = basesave;
	end
	if length(baseabate)==0
		baseabatement = ones(size(p.t)) * p.miu0;
	elseif length(baseabate)==1
		baseabatement = ones(size(p.t)) * baseabate;
	else
		baseabatement = baseabate;
	end
	
	Nscenarios = length(p.T2xCO2);
	
	basesavings=basesavings(:)';baseabatement=baseabatement(:)';
	[EW,Base] = DICEEconomicModel(p,basesavings,baseabatement,[],incrementalemissions);
	if length(startabate)==0
		startabate=baseabatement;
	end

	Base.savings = basesavings;
	%Base.abatement = baseabatement;

%	Base.Population = p.L;
	Base.times = p.t;

	Base.p = p;
	Base.SCC = SCC(p,basesavings,baseabatement,2);
	Base.exitflag = 1;

	for i=1:length(coarseners)
		if nargout>1	
			Aeq=zeros(length(p.t));
			if coarseners(i)>1
				lockedslots = ones(size(p.t));
				lockedslots(1:coarseners(i):end)=0;
				if length(FixedAbatement)>1
					lockedslots(1:length(FixedAbatement))=0;
				end
				lockedslots(end-1)=0; lockedslots(end)=0;
				subunlocked=find(~lockedslots);
				index=1;
				for j=1:length(lockedslots)
					if lockedslots(j)
						pos = j-prevunlocked;
						if nextunlocked<=length(p.t)
							wts = [(1-pos/(nextunlocked-prevunlocked)) (pos/(nextunlocked-prevunlocked))];
							Aeq(j,[prevunlocked nextunlocked]) = [wts];
						else
							Aeq(j,[prevunlocked]) = [1];
						end			
					else
						prevunlocked = subunlocked(index);
						if index<length(subunlocked)
							nextunlocked = subunlocked(index+1);
						else
							nextunlocked = prevunlocked+coarseners(i);
						end
						index=index+1;
						Aeq(j,j) = 1;
					end
				end
			else
				lockedslots=zeros(size(p.t));
				Aeq=eye(length(p.t));
				subunlocked=find(~lockedslots);
			end



			% first optimize only abatement
			WelfAbate = @(x) Welfare([basesavings(:) ; Aeq(:,subunlocked)*x]',p,incrementalemissions);
			
			monotonicityA = sparse(zeros(length(subunlocked)));
			if doMonotonicity
				for i=1:length(subunlocked)-1
					monotonicityA(i,i:i+1) = [1 -1];
				end
			end
			monotonicityB = sparse(zeros(length(subunlocked),1));
			
			startingPoint = [startabate(:)' ]';
	%		Lwr = [[p.miu_2005 zeros(1,timesteps-1)] ]';
	%		Upr = [[p.miu_2005 (ones(1,timesteps-1) * p.limmiu)] ]';
			Lwr = [[ zeros(1,timesteps)] ]';
			if length(p.limmiu)==1
				Upr = [[ (ones(1,timesteps) * p.limmiu)] ]';
			else
				Upr = p.limmiu(:);
			end
			Lwr(1:length(FixedAbatement))=FixedAbatement;
			Upr(1:length(FixedAbatement))=FixedAbatement;

			if usenlopt
			% experimental code
				fitoptions.min_objective = @(x) WelfAbate(x');
				fitoptions.verbose = 1;
				%fitoptions.xtol_rel = 1e-3;
				%fitoptions.ftol_rel = 1e-8;
				fitoptions.xtol_abs = [ones(1,length(startingPoint))*1e-3];
				fitoptions.maxeval = MaxFunEvals;
				fitoptions.algorithm = NLOPT_LN_BOBYQA;
				fitoptions.lower_bounds = Lwr(subunlocked);
				fitoptions.upper_bounds = Upr(subunlocked);
				[xopt,fopt,retcode] = nlopt_optimize(fitoptions,startingPoint);
				xopt = Aeq(:,subunlocked)*xopt;
				xopt = xopt(:)';
			else		
			
	%			fitoptions=optimset('Display',dispmode,'TolX',TolX,'TolFun',abs(1e-8*Base.ExpectedWelfare),'MaxFunEval',MaxFunEvals,'Algorithm',algorithm);
				fitoptions=optimset('Display',dispmode,'MaxFunEval',MaxFunEvals,'Algorithm',algorithm,'UseParallel',useparallel,'TolX',TolX,'TolFun',TolFun);
	
				[optm1.coeffs,optm1.fval,optm1.exitflag,optm1.output,optm1.lambda,optm1.grad] = fmincon(WelfAbate,startingPoint(subunlocked),full(monotonicityA),full(monotonicityB),[],[],Lwr(subunlocked),Upr(subunlocked),[],fitoptions);
				xopt = Aeq(:,subunlocked)*optm1.coeffs;
				xopt = xopt(:)';
			end
	
		[EW,OptimAbatement]	= DICEEconomicModel(p,basesavings(:)',xopt,[],incrementalemissions);
				%OptimAbatement.abatement=optm1.coeffs;
			OptimAbatement.savings=basesavings;	
		
			OptimAbatement.times = p.t;
	
			OptimAbatement.p = p;
			OptimAbatement.SCC = SCC(p,OptimAbatement.savings,OptimAbatement.abatement,2);
			OptimAbatement.exitflag = optm1.exitflag;
			
			baseabatement = xopt;
	
		end

		if nargout > 2	
			
			Aeq2 = [Aeq Aeq*0 ; Aeq*0 Aeq];
			subunlocked2 = [subunlocked subunlocked+length(p.t)];
			% then optimize savings and abatement
			Welfare2 = @(x) Welfare([Aeq2(:,subunlocked2) * x(:)]',p,incrementalemissions);
		
			monotonicityA2 = sparse(zeros(length(subunlocked2)));
			monotonicityA2(length(subunlocked)+1:end,length(subunlocked)+1:end)=monotonicityA;
			monotonicityB2 = sparse(zeros(length(subunlocked2),1));
		
                        basesavings=basesavings(:)';
			startingPoint = [basesavings(subunlocked) optm1.coeffs(:)' ]';
			Lwr = [[basesavings(1) ones(1,timesteps-1)*.01] [p.miu0 zeros(1,timesteps-1)] ]';
			if length(p.limmiu)==1
				Upr = [[basesavings(1) ones(1,timesteps-1)*.9] [p.miu0 (ones(1,timesteps-1) * p.limmiu)] ]';
			else
				Upr = [[basesavings(1) ones(1,timesteps-1)*.9] [p.miu0 p.limmiu(2:end)] ]';
			end
%			fitoptions=optimset('Display',dispmode,'TolX',1e-4,'TolFun',abs(1e-8*Base.ExpectedWelfare),'MaxFunEval',MaxFunEvals,'Algorithm',algorithm);
%			fitoptions=optimset('Display',dispmode,'TolX',1e-6,'TolFun',1e-8,'MaxFunEval',MaxFunEvals,'Algorithm',algorithm);

[optm2.coeffs,optm2.fval,optm2.exitflag,optm2.output,optm2.lambda,optm2.grad] = fmincon(Welfare2,startingPoint,full(monotonicityA2),full(monotonicityB2),[],[],Lwr(subunlocked2),Upr(subunlocked2),[],fitoptions);
				optm2.coeffs = Aeq2(:,subunlocked2)*optm2.coeffs;
				optm2.coeffs = optm2.coeffs(:)';
	
		[EW,FullOptim] = DICEEconomicModel(p,optm2.coeffs(1:timesteps),optm2.coeffs(timesteps+1:end),[],incrementalemissions);
		
			%FullOptim.abatement = [optm2.coeffs(timesteps+1:end)];
			FullOptim.savings = [optm2.coeffs(1:timesteps)];
		
%			FullOptim.Population = p.L;
			FullOptim.times = p.t;
			
			FullOptim.p = p;
		
			FullOptim.SCC = SCC(p,FullOptim.savings,FullOptim.abatement,2);
			FullOptim.exitflag = optm2.exitflag;
			basesavings = FullOptim.savings;
			baseabatement = optm2.coeffs;
			
		end
	
	end

end	

function wlf=Welfare(x,p,incrementalemissions)

	x=x(:)';
	lentime = length(p.t);
	wlf = (-1 * DICEEconomicModel(p,x(1:lentime),x(lentime+1:end),[],incrementalemissions));

end