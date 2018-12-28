function [F] = faceFluxes(casedef)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
dom = casedef.dom;
u = casedef.U.data(1,:);
v = casedef.U.data(2,:);
F = zeros(dom.nC,1);
for i= 1:dom.nF
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    Af = dom.fArea(i);
    lambda = getLambda(dom,i);
    n = dom.fNormal(:,i);
    Uf = [lambda*u(firstCell) + (1-lambda)*u(secondCell); ...
        lambda*v(firstCell) + (1-lambda)*v(secondCell)];
    outwardFlux = Af*Uf'*n;
    if i > dom.nIf
        F(firstCell) = F(firstCell);
    end
    F(firstCell) = F(firstCell) + outwardFlux;
    if secondCell <= dom.nPc
        F(secondCell) = F(secondCell) - outwardFlux;
    end
end

end

