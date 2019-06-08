function [Fo1, Eo, Ed, LatentHeat] = find_e( q )
%FIND_E Work out the Eo, Ed values (and subsequently internal energy and
%   latentheat for a given q value

%ed+eo=y
edpeo = -2*(1+1/sqrt(q));

%2cosh(theta) =sqrt(q)
theta = acosh(sqrt(q)/2);

% LatentHeat =ed-eo=x
LatentHeat = 2*(1+1/sqrt(q))*tanh(theta/2)*tanh(theta)^2;
z=1/sqrt(q);
Fo1 = (1+z)+LatentHeat/2;

Ed = (LatentHeat + edpeo) / 2;
Eo = edpeo - Ed;

end

