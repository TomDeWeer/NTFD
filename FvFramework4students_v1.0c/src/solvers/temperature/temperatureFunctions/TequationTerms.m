function [anb,firstConvTerm,secondConvTerm] ...
    = tEquationTerms( casedef, faceIndex )
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
anb = kappa*Af/Lxi; % must be added (not subtracted) to diagonal terms, because diffusive flux is positive from 1 to 2
% Calculating convection terms of the equations
n = dom.fNormal(:,faceIndex);
u = (lambda*U.data(:,firstNbC) + (1-lambda)*U.data(:, secondNbC)); % works if U is defined on the cells
% u = U.data(:,i); % works if U is defined on the faces
unf = dot(u,n); % n's direction is from firstNbC to second -> convective flux to 1 is positive
% firstConvTerm*phi_1 is the convective contribution from the given face of
% phi at the first cell. It is therefore a diagonal term.
% -secondConvTerm*phi_2 is the contribution of the second cell to the
% equation of the first cell.
% Why the minus sign? dphi_conv_face->cell1 = unf*phi_f and phi_f = lambda*phi_1+(1-lambda)*phi_2
% Since secondConvTerm is later subtracted, the extra minus here is
% necessary
firstConvTerm = lambda*unf*Af;      % the +- signs are
secondConvTerm = -(1-lambda)*unf*Af; % determined by n
        
end

