function ple = double_gaussian( T, E, q, L, scale, Eo, Ed, Co, Cd)
%DOUBLE_GAUSSIAN Double gaussian approximation to density of states
% for single T value, over a range of E values at a given q and lattice
% size. Only q=5,7,10 work. DoS File is used simply to quickly work out
% E+ and E-. Should probably come up with an alternate approach

d = 2;
Tc = 1/log(1+sqrt(q));
kb = 1;
if nargin < 8
    % fprintf('Calculating Cd, Co\n');
    [Cd, Co] = get_quantities(q);
else
    fprintf('Was passed Cd, Co: %f, %f\n', Cd, Co);
end


dt = T-Tc;
adenom = 2*kb*T*Tc;
anumer = dt*(Ed-Eo)*(L^d);
ad = exp(anumer/adenom);
ao = q*exp(-anumer/adenom);

tkb2 = 2*kb*T*T;
ple = scale.*((ad/sqrt(Cd)).*exp(-(((E-Ed-Cd*dt).^2).*(L^d))./(tkb2*Cd))+(ao/sqrt(Co)).*exp(-(((E-Eo-Co*dt).^2).*(L^d))./(tkb2*Co)));

end

function [Cd, Co] = get_quantities(q)
    
    if q==5
        Cd = 2889.2;
        Co = 2886.3;
    elseif q==7
        Cd = 68.7382;
        Co = 68.5132;
    elseif q==10
        %Arisue's values
        Cd = 18.385432;
        Co = 17.937802;
        % Brigg's values
        %Cd = 33; %pm 3
        %Co = 31.8; %pm 2.8
        % Test
        %Cd = 3.4;
        %Co = 3.1;
    else
        fprintf('q must be 5,7 or 10, was: %d', q);
    end
end
