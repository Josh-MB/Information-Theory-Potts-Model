function H = entropy(p)

% Calculate H(X) where p is the probability density of the (possibly
% multi-dimensional) random variable X.

H = -sum(xlogx(p(:)));

function h = xlogx(p)

h = p.*log(p);
h(p < eps) = 0;
h(p < -eps) = NaN;
