function [viscTerm,firstConvTerm,secondConvTerm,pressureForce] ...
    = mgEquationTerms( casedef, faceIndex )
%EQUATIONTERMS - Gives the equation terms for a given face.

dom = casedef.dom;
u = casedef.U.data(1,:);
v = casedef.U.data(2,:);
p = casedef.P.data;

% Determining the indices of the neighbouring cells
[firstNbC,secondNbC] = getCells(dom,faceIndex);
lambda = getLambda(dom,faceIndex);
% Calculating diffusion terms of the equations
Af = dom.fArea(faceIndex);
Lxi = dom.fXiMag(faceIndex);
nu = casedef.material.nu;
viscTerm = nu*Af/Lxi; 
% Calculating convection terms of the equations
n = dom.fNormal(:,faceIndex);
uface = [lambda*u(firstNbC) + (1-lambda)*u(secondNbC); ...
   lambda*v(firstNbC) + (1-lambda)*v(secondNbC)]; % works if U is defined on the cells
% u = U.data(:,i); % works if U is defined on the faces
unf = dot(uface,n); % n's direction is from firstNbC to second -> convective flux to 1 is positive
% Computing pressure contribution
pface = lambda*p(firstNbC) + (1-lambda)*p(secondNbC);
rho = casedef.material.rho;
pressureForce = -Af*pface*n/rho;

firstConvTerm = lambda*unf*Af;      % the +- signs are
secondConvTerm = -(1-lambda)*unf*Af; % determined by n
        
end

