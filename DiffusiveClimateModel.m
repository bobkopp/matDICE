function [Taccum]=DiffusiveClimateModel(Forcing,Tinits,climsens,ocean_diffusivity,sbconst,Cp_mixed,z_mixed,z_deep,Ndeepboxes,dt,F2xCO2)

% [T]=DiffusiveClimateModel(Forcing,Tinits,climsens,sbconst,Cp_mixed,ocean_diffusivity,z_mixed,z_deep,Ndeepboxes,dt,F2xCO2)
%
% Last updated by Robert E. Kopp rkopp-at-alumni.caltech.edu, 17 April 2012


	defval('F2xCO2',3.7127);
	defval('climsens',3);
	defval('sbconst',5.7e-8 * (pi*1e7));
	defval('focean',0.7);
	defval('z_mixed',100);
	defval('Cp_mixed',focean * 3985 * 1030 * z_mixed);
	defval('z_deep',3790-z_mixed);
	defval('Ndeepboxes',2);
	defval('dt',1);
	defval('ocean_diffusivity',1e-3 * (pi*1e7));

	defval('Tinits',zeros(1,Ndeepboxes+1));
	
	lam=F2xCO2./climsens;

	zs = [z_mixed repmat(z_deep/Ndeepboxes,1,Ndeepboxes)];
	dists = .5*(zs(1:end-1)+zs(2:end));
	mixfluxcoeffs = bsxfun(@rdivide,ocean_diffusivity,dists);
	zmin = min([zs(1:end-1);zs(2:end)],[],1);
	
	Nsow = size(Forcing,1);
	
	if size(Tinits,1) ~= Nsow
		Tinits = repmat(Tinits,Nsow,1);
	end
	
	% in Ts: rows are States of the World, columns are boxes, 3rd dim is time
	
	Taccum(:,:,1) = Tinits;

	for j=1:size(Forcing,2)
		clear Ts;
		Ts(:,:,1) = Taccum(:,:,j);
		Teq = bsxfun(@rdivide,Forcing(:,j),lam);
		for i=2:(dt+1)

			Ts(:,:,i) = Ts(:,:,i-1);
			dTs = Ts(:,1:end-1,i)-Ts(:,2:end,i);

			forcedwarm = sbconst./Cp_mixed * (4*255^3) * (Teq - Ts(:,1,i-1));
			Ts(:,1,i) = Ts(:,1,i) + forcedwarm;
			dTs2 = Ts(:,1:end-1,i)-Ts(:,2:end,i);
			mixfluxes = bsxfun(@times,mixfluxcoeffs,.5*(dTs+dTs2));

			mixfluxes=sign(mixfluxes).*min(abs(bsxfun(@times,zmin,dTs2)),abs(mixfluxes));
			Ts(:,1:end-1,i) = Ts(:,1:end-1,i) + bsxfun(@rdivide,-mixfluxes,zs(1:end-1));
			Ts(:,2:end,i) = Ts(:,2:end,i) + bsxfun(@rdivide,mixfluxes,zs(2:end));
			dTs3 = Ts(:,1:end-1,i)-Ts(:,2:end,i);
	


%			if sum((dTs2./dTs)<0)>0
%				keyboard
%			end
%				
%			Ts(:,:,i) = min(Ts(:,:,i),max(Ts(:,:,i-1)))
%
%			signflippers = ( (Ts(:,1,i)-Ts(:,2,i))./(Ts(:,1,i-1)-Ts(:,2,i-1)+1e-5) ) < 0;
%			if sum(signflippers)>0
%				keyboard
%			end
%			Ts(:,1,i) = ~signflippers .* Ts(:,1,i) + signflippers .* Ts(:,2,i);
%
%			signflippers = ( (Ts(:,end,i)-Ts(:,end-1,i))./(Ts(:,end,i-1)-Ts(:,end-1,i-1)+1e-5) ) < 0;
%			Ts(:,1,i) = ~signflippers .* Ts(:,end,i) + signflippers .* Ts(:,end-1,i);
%
%			Ts(:,1,i) = Ts(:,1,i) + forcedwarm;
		end
		Taccum(:,:,j+1) = Ts(:,:,end);
	end

	% now reshape so that rows are time, columns are boxes, 3rd dim is SOW
	Taccum = permute(Taccum,[3 2 1]);

end
