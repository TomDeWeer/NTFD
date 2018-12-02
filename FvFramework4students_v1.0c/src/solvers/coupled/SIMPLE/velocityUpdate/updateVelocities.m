function [Uupdated] = updateVelocities(casedef,Pcorr,uP,vP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dom = casedef.dom;

% Create field objects
U = Field(dom.allCells,1);	% Velocity [m/s] (2D vector); empty field
reset(U,[0; 0]); 
set(U,casedef.U.data);      % Set to given initial guess

u = U.data(1,:);
v = U.data(2,:);

for i= 1:nIf+nBf
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    % Computing pressure contribution
    pCorrFace = lambda*Pcorr(firstNbC) + (1-lambda)*Pcorr(secondNbC);
    rho = casedef.material.rho;
    pressureCorrForce = -Af*pCorrFace*n/rho;
    % First cell
    u(firstCell) = u(firstCell) + pressureCorrForce/uP(firstCell);
    v(firstCell) = v(firstCell) + pressureCorrForce/vP(firstCell);
    % Second cell
    % Checking whether it's a physical cell or a ghost cell
    if secondCell <= dom.nPc % Physical cell
        u(secondCell) = u(secondCell) - pressureCorrForce/uP(secondCell);
        v(secondCell) = v(secondCell) - pressureCorrForce/vP(secondCell);
    else % If it's a ghost cell
        [uGhost, vGhost] = UcorrGhostTerms(casedef, i, u(firstCell), v(firstCell));
        u(secondCell) = uGhost;
        v(secondCell) = vGhost;
    end
end

Uupdated = [u, v];

end

