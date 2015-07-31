function R = icdftri(seeds,exptri)

% R = icdftri(seeds,[a b c])
%
% Takes seeds (between 0 and 1) and inverts the CDF from a triangular distribution
% with minimum a, mode b, and maximum c.
%
% Last updated by Robert E. Kopp rkopp-at-alumni.caltech.edu, 22 May 2011
% Copyright (C) 2011 by Robert E. Kopp; distributed under GNU GPL v3

defval('exptri',[1 1.5 3.5]);
a=exptri(1); b=exptri(2); c=exptri(3);
trs = ( b - a ) ./ ( c - a); 

R = (seeds<=trs).*(a + sqrt( (c-a) .* (b-a) .* seeds)) + (seeds>trs) .* ( c - sqrt ( (1 - seeds) .* (c - a) .* ( c - b ) ) );