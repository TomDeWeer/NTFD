function [Pcorr] = pressureCorrSolver(casedef, uP, vP)
%pressureCorrSolver Summary of this function goes here
%   Detailed explanation goes here

dom = casedef.dom;

Pcorr = Field(dom.allCells,0);	% Pressure [Pa] (scalar); empty field
set(Pcorr,zeros(1,dom.nC));

A = sparse(double(dom.nC),double(dom.nC)) ; % contains pressure correction equations(1 for every pressure, even ghostcells)

F = faceFluxes(casedef);      % Without Rie-Chow
% F = faceFluxesRC(casedef, uP, vP);  % With Rie-Chow

% Creating pressure corrections equations in internal cells:
for i= 1:dom.nIf+dom.nBf
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    lambda = getLambda(dom,i);
    n = dom.fNormal(:,i);
    Af = dom.fArea(i);
    %%%%% TODO: TOM SNAPT DIT NIET (Koen eigenlijk ook niet) %%%%%
%     theta = atan2(n(2),n(1));
%     if secondCell <= dom.nPc
%         af = sqrt(cos(theta)^2*(lambda*uP(firstCell) + (1-lambda)*uP(secondCell))^2 ...
%             + sin(theta)^2*(lambda*vP(firstCell) + (1-lambda)*vP(secondCell))^2);
%     else %nu heb je geen tweede vergelijking
%         af = sqrt(cos(theta)^2*uP(firstCell)^2 + sin(theta)^2*vP(firstCell)^2);
%     end
    if secondCell <= dom.nPc
        af = lambda*uP(firstCell) + (1-lambda)*uP(secondCell);
    else %nu heb je geen tweede vergelijking
        af = uP(firstCell);
    end
    df = Af^2/af;
    % op randfaces heb je niet langs beide kanten een momentumvgl -> pak
    % alleen die van internal cell
    A(firstCell,firstCell) = A(firstCell,firstCell) + df; % equation for firstCell
    A(firstCell,secondCell) = A(firstCell,secondCell) - df; % equation for firstCell
    if secondCell <= dom.nPc
        A(secondCell,secondCell) = A(secondCell, secondCell) + df; % equation for secondCell
        A(secondCell,firstCell) = A(secondCell, firstCell) - df; % equation for secondCell
    end
end


% add pressure correction boundary conditions
for faceIndex= dom.nIf+1:dom.nF
    [physicalCell,ghostCell] = getCells(dom,faceIndex);
    current_ghost_p = casedef.P.data(ghostCell);
    current_physical_p = casedef.P.data(physicalCell);
    % Checking which boundary the face belongs to
    for randID = 1:length(casedef.BC)
        range = dom.getzone(casedef.BC{randID}.zoneID).range;
        if faceIndex >= range(1) && faceIndex <= range(end)
            id =  randID;
            break
        end
    end
    % Checking which BC applies at that boundary
    BC = casedef.BC{id}.pressureKind;
    switch BC
        case 'Dirichlet' % correction term must make sure that dirichlet condition is satisfied
            % Determining lambda using the anonymous function
            lambda = getLambda(dom,faceIndex);
            p_hat = casedef.BC{id}.data.pressure;
            if isa(p_hat, 'function_handle')
                pos = dom.fCoord(:,faceIndex);
                p_described = p_hat(pos(1),pos(2));
            else
                p_described = p_hat;
            end
            A(ghostCell,ghostCell) = lambda;
            A(ghostCell,physicalCell) = 1-lambda;
            F(ghostCell) = -(p_described-lambda*current_physical_p-(1-lambda)*current_ghost_p);
        case 'Neumann'
            ksi = dom.fXiMag(faceIndex);
            dp = casedef.BC{id}.data.pressure;
            if isa(dp, 'function_handle')
                pos = dom.fCoord(:,faceIndex);
                dp_described = dp(pos(1),pos(2));
            else
                dp_described = dp;
            end
            A(ghostCell,ghostCell) = (1/ksi);
            A(ghostCell,physicalCell) = -(1/ksi);
            F(ghostCell) = -(dp_described - (current_ghost_p-current_physical_p)/ksi);
        otherwise
            disp('BC not found');
    end
end

Pcorr = A\-F;
% PcorrField = Field(dom.allCells,0);
% set(PcorrField, Pcorr')
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("Pcorr"); colorbar();
% fvmplotfield(PcorrField,scale,lw);



end

