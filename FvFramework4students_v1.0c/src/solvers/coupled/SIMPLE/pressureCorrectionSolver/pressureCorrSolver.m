function [Pcorr] = pressureCorrSolver(casedef, uP, vP)
%pressureCorrSolver Summary of this function goes here
%   Detailed explanation goes here

dom = casedef.dom;

Pcorr = Field(dom.allCells,0);	% Pressure [Pa] (scalar); empty field
set(Pcorr,zeros(1,dom.nC));

A = sparse(double(dom.nC),double(dom.nC)) ; % contains pressure correction equations(1 for every pressure, even ghostcells)

F = faceFluxes(casedef);
F = F + RCcorrectionTerms(casedef, uP, vP);  % Rhie-Chow correction

% Creating pressure corrections equations in internal cells:
afGrad = zeros(1,dom.nC);
for i= 1:dom.nIf+dom.nBf
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    lambda = getLambda(dom,i);
    Af = dom.fArea(i);
    ksi = dom.fXiMag(i);
    if secondCell <= dom.nPc
        af = lambda*uP(firstCell) + (1-lambda)*uP(secondCell);
    else %nu heb je geen tweede vergelijking
        af = uP(firstCell);
    end
    rho = casedef.material.rho;
    df = Af^2/(af*rho);
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
            A(ghostCell,ghostCell) = 1-lambda;
            A(ghostCell,physicalCell) = lambda;
            F(ghostCell) = -(p_described-lambda*current_physical_p-(1-lambda)*current_ghost_p);
        case 'Neumann'
            % Check normalization
            if casedef.BC{id}.isNormalized && faceIndex == range(1)
                % apply zero dirichlet BC
                p_described = 0;    % Standard pressure in literature is 0Pa
                A(ghostCell,ghostCell) = 1-lambda;
                A(ghostCell,physicalCell) = lambda;
                F(ghostCell) = -(p_described-lambda*current_physical_p-(1-lambda)*current_ghost_p);
            else
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
            end
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
% 
% FluxField = Field(dom.allCells,0);
% set(FluxField, F')
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("Fluxes"); colorbar();
% fvmplotfield(FluxField,scale,lw);



end

