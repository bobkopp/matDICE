function [TCR,Tatm,Tocean]=DICETransientClimateResponse(c1,c3,c4,lam,FCO22x)

% TCR=DICETransientClimateResponse
%
% Last updated by Robert E. Kopp rkopp-at-alumni.caltech.edu, 16 April 2012

	defp = DICEParameters;
	defval('FCO22x',defp.FCO22x);
	defval('c1',defp.c1);
	defval('c3',defp.c3);
	defval('c4',defp.c4);
	defval('lam',defp.lam);

	calcnexttemperature = @(Atm,Ocean,Forcing) [Atm+c1*(Forcing-lam.*Atm-c3*(Atm-Ocean)) Ocean+c4*(Atm-Ocean)];
	if nargout == 1
		t=0:70;
	else
		t=0:1000;
	end
	ppmCO2 = 280 * 1.01.^t;
	ppmCO2(find(ppmCO2>560))=560;
	RF = log(ppmCO2/280)/log(2) * FCO22x;
	
	N=max([length(c1) length(c3) length(c4) length(lam)]);
	Tatm(:,1) = zeros(N,1);
	Tocean(:,1) = zeros(N,1);
	for i=2:length(t)
		temps = calcnexttemperature(Tatm(:,i-1),Tocean(:,i-1),RF(i));
		Tatm(:,i) = temps(:,1); Tocean(:,i) = temps(:,2);
	end
	TCR=Tatm(:,71);
	
	
end
		