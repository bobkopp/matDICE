function y = pmean(X,p)

% Last updated by Bob Kopp rkopp-at-alumni.caltech.edu, 9 June 2011

if p < 0
	m = min(abs(X(:)));
else
	m = max(abs(X(:)));
end

X2 = X/m;

q = mean((X2.^p));
y = (q.^(1./p))*m;
