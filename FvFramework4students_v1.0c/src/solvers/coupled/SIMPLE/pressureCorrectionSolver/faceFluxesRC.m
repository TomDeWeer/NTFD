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
    omega1 = dom.cVol(firstCell); 
    omega2 = dom.cVol(secondCell); 
    Af = dom.fArea(faceindex);
    n = dom.fNormal(:,faceindex);
    lambda = getLambda(dom,faceindex);
    RCterm =  n*Af*(lambda*P(firstCell)+(1-lambda)*P(secondCell));
    gradP(:,firstCell) = gradP(:,firstCell) + RCterm/omega1;
    %%%% Iets deftigs voor gradP in ghostcells
    if secondCell <= dom.nPc
        gradP(:,secondCell) = gradP(:,secondCell) - RCterm/omega2;
    else %nu heb je geen tweede vergelijking
        gradP(:,secondCell) = gradP(:,firstCell);
    end
end
% gradPfield = Field(casedef.dom.allCells,1);
% set(gradPfield,gradP);
% casedef.gradP = gradPfield;
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("P"); colorbar();
% fvmplotfield(casedef.P,scale,lw);
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("gradPx"); colorbar();
% fvmplotfield(casedef.gradP,scale,lw,1);
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("gradPy"); colorbar();
% fvmplotfield(casedef.gradP,scale,lw,2);
F = zeros(dom.nC,1);
for i= 1:dom.nIf+dom.nBf
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    lambda = getLambda(dom,i);
    Af = dom.fArea(i);
    omega1 = dom.cVol(firstCell); 
    omega2 = dom.cVol(secondCell); 
    n = dom.fNormal(:,i);
    theta = atan2(n(2),n(1));
    if secondCell <= dom.nPc
        af = lambda*uP(firstCell) + (1-lambda)*uP(secondCell);
    else %nu heb je geen tweede vergelijking
        af = uP(firstCell);
    end
    ksi = dom.fXiMag(faceindex);
    directDP = (P(secondCell)-P(firstCell))/ksi;
    interpolatedDP = (lambda*gradP(:,firstCell) + (1-lambda)*gradP(:,secondCell))'*n;
    RCcorrection = -(directDP - interpolatedDP)/af;
    Uf = [lambda*u(firstCell) + (1-lambda)*u(secondCell); ...
        lambda*v(firstCell) + (1-lambda)*v(secondCell)];
    outwardFlux1 = Af*(Uf'*n+omega1*RCcorrection);
    outwardFlux2 = Af*(Uf'*n+omega2*RCcorrection);
%     if n(1)>0 && i ==1
%         disp(directDP-interpolatedDP)
%     end
    F(firstCell) = F(firstCell) + outwardFlux1;
    F(secondCell) = F(secondCell) - outwardFlux2; % kijk hier een uur na of je gewoon min mag doen
end

end

