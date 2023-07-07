function u = matrixcalculator(matxy,par)
%---------------------------MATERIAL PROPERTIES --------------------------%
%-------------------------------------------------------------------------%
EE = 160000;            % YOUNG MODULUS   [MPa]
ni = 0.3;               % POISSON COEFFICIENT
GG = EE/(2*(1+ni));     % SHEAR MODULUS

%kk = (3-4*ni) ;        % PLAIN STRAIN
kk = (3-ni)/(1+ni);     % PLAIN STRESS
%-------------------------------------------------------------------------%

% Model parameters
L = length(par);
a = par(1:(L-4)/2);    % Coefficient for the first mode expansion
b = par((L-4)/2+1:L);  % Coefficient for the second mode expansion
sX = par(L-3);         % Translation along the x-axis
sY = par(L-2);         % Translation along the y-axis
sZ = par(L-1);         % Translation along the z-axis
alpha = par(L);        % Rotation around the z-axis

% Number of terms of the Williams' expansion
n = length(a);

% Initialization with the Z offset
uX = sZ;
uY = sZ;

% X and Y offset 
X = matxy(:,1) + sX;
Y = matxy(:,2) + sY;

% Conversion to polar coordinates
r = sqrt(X.^2 + Y.^2);
th = atan2(Y,X) - alpha ; % If alpha > 0: counterclockwise rotation

% X and Y displacement by William's expansion
for i = 1:n
    uX = uX + r.^(i/2)/(2*GG).*a(i).*((kk+i/2+(-1)^i).*cos(i*th/2)-i/2.*cos(th*(i-4)/2)) ...
            - r.^(i/2)/(2*GG).*b(i).*((kk+i/2-(-1)^i).*sin(i*th/2)-i/2.*cos(th*(i-4)/2));
        
    uY = uY + r.^(i/2)/(2*GG).*a(i).*((kk-i/2-(-1)^i).*sin(i*th/2)+i/2.*sin(th*(i-4)/2)) ...
            + r.^(i/2)/(2*GG).*b(i).*((kk-i/2+(-1)^i).*cos(i*th/2)+i/2.*cos(th*(i-4)/2));
end

u = uY; % Replace this line with "u = uX" if X displacement is needed

end