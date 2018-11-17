function [CE] = NavierStokes(casedef, x)
% Computes the residuals of the discretized solution to the Navier Stokes
% equations, and the Jacobian.

dom = casedef.dom;

% unpack x
[p, u, v] = getPUV(casedef, x);

% Continuity equations: one for every physical cell, contributions for
% every face
CE = zeros(dom.nPc,1);

for i=1:dom.nF
    Af = dom.fArea(i);
    [firstNbC,secondNbC] = getCells(dom,i);
    lambda = getLambda(dom,i); % points from firstNbc to secondNbc
    u1 = [u(firstNbC); v(firstNbC)];
    u2 = [u(secondNbC); v(secondNbC)];
    % equation for first Nbc
    Xi = dom.fXi(:,i); % pointing correctly for firstNbc
    ni = Xi/norm(Xi,2);
    CE(firstNbC) = CE(firstNbC) + Af*lambda*ni'*u1 + Af*(1-lambda)*ni'*u2;
    if secondNbC<=dom.nPc
        CE(secondNbC) = CE(secondNbC) - Af*lambda*ni'*u1 - Af*(1-lambda)*ni'*u2;
    end
end






end
