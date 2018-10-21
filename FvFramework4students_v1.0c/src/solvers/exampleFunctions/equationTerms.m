function [anb,firstConvTerm,secondConvTerm] ...
    = equationTerms( casedef, faceIndex )
%EQUATIONTERMS - Gives the equation terms for a given face.
%Function to determine the terms of the internal convection-diffusion
%equations. 'anb' is the diffusive term. 'firstConvTerm' ('secondConvTerm')
%is the convective term of the first (second) cell. Both should be added to
%their respective diagonal matrix elements, and subtracted from the other's
%off-diagonal row element.

dom = casedef.dom;
U = casedef.U0;

% Determining the indeces of the neighbouring cells
[firstNbC,secondNbC] = getCells(dom,faceIndex);
lambda = getLambda(dom,faceIndex);
% Calculating diffusion terms of the equations
Af = dom.fArea(faceIndex);
Lxi = dom.fXiMag(faceIndex);
kappa = casedef.material.k;
anb = kappa*Af/Lxi;
% Calculating convection terms of the equations
n = dom.fNormal(:,faceIndex);
u = (lambda*U.data(:,firstNbC) + (1-lambda)*U.data(:, secondNbC)); % works if U is defined on the cells
% u = U.data(:,i); % works if U is defined on the faces
unf = dot(u,n); % n's direction is from firstNbC to second
firstConvTerm = lambda*unf*Af;      % the +- signs are
secondConvTerm = -(1-lambda)*unf*Af; % determined by n
        
end

