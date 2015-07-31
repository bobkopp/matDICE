function [spec,seeds]=DamFuncXau(N,uncscalefac,seeds,cses,aa2_gradual_calib,aa3_mean,aa3_logvar,aa2_mean,aa2_logvar,aa_cat_damage_mean,aa_cat_damage_var,aa_cat_width_mean,aa_cat_width_logvar)

% y=DamFuncXau(scalefactor)
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Jul 14 21:53:27 EDT 2012


defp=DICEParameters;

defval('N',100)
defval('uncscalefac',1);
defval('seeds',[]);
defval('cses',[]);
defval('aa2_gradual_calib',.0061);
defval('aa3_mean',2);
defval('aa3_logvar',log(2)/2);
defval('aa2_mean',defp.aa2);
defval('aa2_logvar',log(2));
defval('aa_cat_damage_mean',0.3);
defval('aa_cat_damage_var',.15);
defval('aa_cat_width_mean',0.1);
defval('aa_cat_width_logvar',log(2));

if length(seeds)<N
% generate some minimally correlated "seeds" (vector of randon numbers from 0 to 1) for use in Latin Hypercube sampling of uncertain parameters

	for j=1:20 % number of samples to draw and check for correlation
		seq(:,1,j)=1:N;
		for i=2:10
			seq(:,i,j)=randperm(N);
		end
		cor = corrcoef(seq(:,:,j)).^2;
		cors(j) = sum(cor(:));
	end
	[m,mi]=min(cors); % find the least correlated set of seeds
	seeds = seq(:,:,mi)/N - rand(size(seq,1),size(seq,2))/N; % take that one and normalize so that it maps to [0,1]

end

if length(cses)==1
	cses = ones(N,1) * cses;
elseif length(cses)<1
	cses = icdfRoeBaker(seeds(:,1));
end
	

% Now we start defining different damage functions [this is a fairly long section]
% for each damage function, we store the key parameters we'll later pass on to DICEParameters in the cell array spec_dam{i,j}, where
% i is the index of the damage function, and j refers to the different calibrations of the damage function defined above.

	% Xbu with variable cat width

	icdfnorm = @(x,mu,v) (-sqrt(2).*erfcinv(2*x).*sqrt(v) + mu);

	aa3_logvar = uncscalefac*aa3_logvar;
	aa3_mean_lognormal = log(aa3_mean) - aa3_logvar/2;
	aa3 = exp(icdfnorm(seeds(:,2),aa3_mean_lognormal,aa3_logvar));

	aa2_logvar = uncscalefac*aa2_logvar;
	aa2_mean_lognormal = log(aa2_mean) - aa2_logvar/2;
	aa2 = exp(icdfnorm(seeds(:,3),aa2_mean_lognormal,aa2_logvar));
	aa2_gradual = (1./(1-aa2_gradual_calib*(aa2/defp.aa2)) - 1) ./ 2.5.^2;
	
	aa2 = aa2 .* 2.5.^(defp.aa3-aa3);
	aa2_gradual = aa2_gradual .* 2.5.^(defp.aa3-aa3);
	aa2_gradual=min(aa2,aa2_gradual);

	aa_cat_damage = aa_cat_damage_mean + uncscalefac*aa_cat_damage_var*icdftri(seeds(:,4),[-1 0 1]); 
	aa2_cat_damage = aa_cat_damage./(1-aa_cat_damage); % coefficient for the denominator
	
	aa_cat_width_logvar = uncscalefac*aa_cat_width_logvar;
	aa_cat_width_lognormal = log(aa_cat_width_mean) - aa_cat_width_logvar/2;
	aa_cat_width = exp(icdfnorm(seeds(:,5),aa_cat_width_lognormal,aa_cat_width_logvar));
	
	% ICDF for temperature threshold for catastrophe
	
	cat_prob = @(Temperature) (1./(1+aa2.*abs(Temperature).^aa3) - 1./(1+aa2_gradual.*abs(Temperature).^aa3))./(-1./(1+aa2_gradual.*abs(Temperature).^aa3) + 1./(1+aa2_gradual.*abs(Temperature).^aa3 + aa2_cat_damage)); 
	
	testT = 0:.1:50;
	clear testT_cat_prob;
	for i=1:length(testT)
		testT_cat_prob(:,i)=cat_prob(testT(i));
	end

	clear cat_thresholds;
	for i=1:N
		[u,ui]=unique(testT_cat_prob(i,:));
		if length(u)>1
			cat_thresholds(i) = interp1(testT_cat_prob(i,ui),testT(ui),seeds(i,6),'linear','extrap');
		else
			cat_thresholds(i)=1000;
		end
	end
	
	cat_thresholds=cat_thresholds(:);
	
	spec = {'T2xCO2',cses,'aa3',aa3,'aa2_gradual',aa2_gradual,'aa2_cat_damage',aa2_cat_damage,'calc_aa2_cat_damage_multiplier','@(Temperature) max(1, (1./(1+p.aa2.*abs(Temperature).^p.aa3) - 1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3))./(-1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3) + 1./(1+p.aa2_gradual.*abs(Temperature).^p.aa3 + p.aa2_cat_damage)))','aa_cat_threshold',cat_thresholds,'aa_cat_width',aa_cat_width,'calcdamages','1-1./(1 + p.aa2_gradual.*p.calc_aa2_cat_damage_multiplier(Temperature).*abs(Temperature).^p.aa3 + (1./(1 + exp(-(abs(Temperature)-(p.aa_cat_threshold+5*p.aa_cat_width))./p.aa_cat_width) )).*p.aa2_cat_damage)'};
	