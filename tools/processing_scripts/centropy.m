function H = centropy(p,x)

% Calculate H(X|Y) = H(X,Y)-H(Y) where p is the joint probability density of the
% (possibly multi-dimensional) random variables X,Y and the vector of indices x
% specifies which dimensions of p reference the X variable.

H = entropy(p) - entropy(multisum(p,x)); % second term is entropy of marginal Y distribution

function s = multisum(p,x)

x = x(:);
s = p;
for i=1:length(x)
    s = sum(s,x(i));
end 
s = squeeze(s);
