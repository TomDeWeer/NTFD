function [F] = faceFluxesRC(casedef, uP, vP)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
dom = casedef.dom;
P = casedef.P.data;
u = casedef.U.data(1,:);
v = casedef.U.data(2,:);
% calculating pressure gradient in every cell using FV
gradP = zeros(2,dom.nC);
for faceindex=1:dom.nF
    [firstCell,secondCell] = getCells(dom,faceindex);
    omega = dom.fVolume(faceindex); % of zoiets
    Af = dom.fArea(faceindex);
    n = dom.fNormal(:,faceindex);
    lambda = getLambda(dom,faceindex);
    RCterm =  n*Af*(lambda*P(firstCell)+(1-lambda)*P(secondCell))/omega;
    gradP(:,firstCell) = gradP(:,firstCell) + RCterm;
    %%%% Iets deftigs voor gradP in ghostcells
    if secondCell <= dom.nPc
        gradP(:,secondCell) = gradP(:,secondCell) - RCterm;
    else %nu heb je geen tweede vergelijking
        gradP(:,secondCell) = -gradP(:,firstCell);
    end
end

F = zeros(dom.nC,1);
for i= 1:dom.nIf+dom.nBf
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    Af = dom.fArea(i);
    n = dom.fNormal(:,i);
    theta = atan2(n(2),n(1));
    if secondCell <= dom.nPc
        af = sqrt(cos(theta)^2*(lambda*uP(firstCell) + (1-lambda)*uP(secondCell))^2 ...
            + sin(theta)^2*(lambda*vP(firstCell) + (1-lambda)*vP(secondCell))^2);
    else %nu heb je geen tweede vergelijking
        af = sqrt(cos(theta)^2*uP(firstCell)^2 + sin(theta)^2*vP(firstCell)^2);
    end
    lambda = getLambda(dom,i);
    Uf = [lambda*u(firstCell) + (1-lambda)*u(secondCell); ...
        lambda*v(firstCell) + (1-lambda)*v(secondCell)];
    outwardFlux = Af*Uf'*n;
    ksi = dom.fXiMag(faceindex);
    directDP = (P(secondCell)-P(firstCell))/ksi;
    interpolatedDP = lambda*gradP(firstCell) + (1-lambda)*gradP(secondCell);
    RCcorrection = -cellVolume*(directDP - interpolatedDP)/af;
    F(firstCell) = F(firstCell) + outwardFlux + RCcorrection;
    F(secondCell) = F(secondCell) - outwardFlux - RCcorrection; % kijk hier een uur na of je gewoon min mag doen
end

end

