function [F] = RCcorrectionTerms(casedef, uP, vP)
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
        
%         ksi = dom.fXiMag(faceindex);
%         gradP(:,secondCell) = n*(P(secondCell)-P(firstCell))/ksi;
        
%     % Checking which boundary the face belongs to
%         for randID = 1:length(casedef.BC)
%             range = dom.getzone(casedef.BC{randID}.zoneID).range;
%             if faceindex >= range(1) && faceindex <= range(end)
%                 id =  randID;
%                 break
%             end
%         end
%         % Checking which BC applies at that boundary
%         BC = casedef.BC{id}.pressureKind;
%         switch BC
%             case 'Dirichlet' % correction term must make sure that dirichlet condition is satisfied
%                 % Determining lambda using the anonymous function
%                 p_hat = casedef.BC{id}.data.pressure;
%                 if isa(p_hat, 'function_handle')
%                     pos = dom.fCoord(:,faceindex);
%                     p_described = p_hat(pos(1),pos(2));
%                 else
%                     p_described = p_hat;
%                 end
%                 for oppositeFace=1:dom.nIf
%                     [cell1,cell2] = getCells(dom,oppositeFace);
%                     extraNormal = dom.fNormal(:,oppositeFace);
%                     if cell1 == firstCell && isequal(n,-extraNormal)
%                         extraCell = cell2;
%                     elseif cell2 == firstCell && isequal(n,-extraNormal)
%                         extraCell = cell1;
%                     else
%                         continue;
%                     end
%                     distance = norm(dom.cCoord(:,extraCell)-dom.cCoord(:,secondCell));
%                     gradP(:,secondCell) = -n*(P(secondCell) - P(extraCell))/distance;
%                     break;
%                 end
% 
%             case 'Neumann'
%                 dp = casedef.BC{id}.data.pressure;
%                 if isa(dp, 'function_handle')
%                     pos = dom.fCoord(:,faceindex);
%                     dp_described = dp(pos(1),pos(2));
%                 else
%                     dp_described = dp;
%                 end
%                 gradP(:,secondCell) = dp_described;
%             otherwise
%                 disp('BC not found');
%         end
        
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
    ksi = dom.fXiMag(i);
    interpolatedDP = (lambda*gradP(:,firstCell) + (1-lambda)*gradP(:,secondCell))'*n;
    directDP = (P(secondCell)-P(firstCell))/ksi;
    if secondCell <= dom.nPc
        af = lambda*uP(firstCell) + (1-lambda)*uP(secondCell);
    else %nu heb je geen tweede vergelijking
        af = uP(firstCell);
    end
    rho = casedef.material.rho;
    af = af*rho;
    RCcorrection = -(directDP - interpolatedDP)/af;
    outwardFlux1 = Af*omega1*RCcorrection;
    outwardFlux2 = Af*omega2*RCcorrection;
    F(firstCell) = F(firstCell) + outwardFlux1;
    F(secondCell) = F(secondCell) - outwardFlux2; % kijk hier een uur na of je gewoon min mag doen
end

end

