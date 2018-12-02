function [Pcorr] = pressureCorrSolver(casedef, uP, vP)
%pressureCorrSolver Summary of this function goes here
%   Detailed explanation goes here

dom = casedef.dom;
u = casedef.U.data(1,:);
v = casedef.U.data(2,:);

Pcorr = Field(dom.allCells,0);	% Pressure [Pa] (scalar); empty field
set(Pcorr,zeros(1,dom.nC));

% Create an equation object for holding a scalar conservation equation
eqnPcorr = ScalarFvEqn2(dom);

F = zeros(dom.nC,1);
A = sparse(double(dom.nC),double(dom.nC)) ; % contains pressure correction equations(1 for every pressure, even ghostcells)

% Creating pressure corrections equations in internal cells:
for i= 1:dom.nIf+dom.nBf
    % Getting terms of the equations
    [firstCell,secondCell] = getCells(dom,i);
    Af = dom.fArea(i);
    lambda = getLambda(dom,i);
    n = dom.fNormal(:,i);
    Uf = [lambda*u(firstCell) + (1-lambda)*u(secondCell); ...
        lambda*v(firstCell) + (1-lambda)*v(secondCell)];
    outwardFlux = Af*Uf'*n;
    F(firstCell) = F(firstCell) + outwardFlux;
    F(secondCell) = F(secondCell) - outwardFlux;
    if secondCell <= dom.nPc
        if n'*[1; 0]>0.5 % face in u richting
            af = (uP(firstCell) + uP(secondCell))/2;
        else
            af = (vP(firstCell) + vP(secondCell))/2;
        end
    else %nu heb je geen tweede vergelijking
        if n'*[1; 0]>0.5 % face in u richting
            af = uP(firstCell); % VRAGEN OF WE MOGEN EXTRAPOLEREN
        else
            af = vP(firstCell);
        end
    end
    % op randfaces heb je niet langs beide kanten een momentumvgl -> pak
    % alleen die van internal cell
    A(firstCell,firstCell) = A(firstCell,firstCell) + af; % equation for firstCell
    A(firstCell,secondCell) = A(firstCell,secondCell) - af; % equation for firstCell
    if secondCell <= dom.nPc
        A(secondCell,secondCell) = A(secondCell, secondCell) + af; % equation for secondCell
        A(secondCell,firstCell) = A(secondCell, firstCell) - af; % equation for secondCell
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
            F(ghostCell) = p_described-lambda*current_physical_p-(1-lambda)*current_ghost_p;
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
            F(ghostCell) = dp_described - (current_ghost_p-current_physical_p)/ksi;
        otherwise
            disp('BC not found');
    end
end

Pcorr = A\F;

end

