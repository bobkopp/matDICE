function [W,U,curu] = EZUtility(C,L,eta,RRA,rho,t)

% [W,U] = EZUtility(C,L,IES,RRA,rho,t)
%
% Epstein-Zin Utility
%
%
% Last updated by Bob Kopp rkopp-at-alumni.caltech.edu, 8 June 2011

defval('L',ones(size(t)));

N = size(C,1);
U = zeros(N,size(C,2));

if eta == 1
	eta = 1.01;
end
if RRA == 1
	RRA = 1.01;
end

alf = 1/(1-eta);
puretime = repmat((1./(1+rho)).^(t-t(1)),N,1);
discfact = puretime.*(repmat(L,N,1));

curu = abs(alf) * C.^(1./alf);
U(:,end) = curu(:,end);
b = (alf*(1-RRA));
for j=(size(C,2)-1):-1:1
	%U(:,j) = curu(:,j) + (discfact(:,j+1)./discfact(:,j)) .*  real((mean ( U(:,j+1).^(alf*(1-RRA)) )).^(1./(alf*(1-RRA))));

%	ubar = mean(U(:,j+1));
%	usigsq = sum((U(:,j+1)-ubar).^2)/N;
%	q1 = log(ubar) + log((1 + b*(b-1)/(ubar^2)*usigsq)/b);
%	q2 = (log(usigsq) + log(.5 * b * (b-1)) + (b-2) * log(ubar) + log (1 + ubar.^2/(b*(b-1).*usigsq)));
%	q = q1;
%	q(find(~isreal(q)+isnan(q)+isinf(q))) = q2;
%	U(:,j) = real(curu(:,j) + (discfact(:,j+1)./discfact(:,j)) .*  exp(q));

	U(:,j) = curu(:,j) + (discfact(:,j+1)./discfact(:,j)) .*  pmean(U(:,j+1),b);
end

W = mean(U(:,1))*L(1);
