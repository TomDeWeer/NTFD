function [ ughostDiag, ughostOffdiag, ubValue, vghostDiag, vghostOffdiag, vbValue ] ...
    = mgGhostTerms( casedef, faceIndex )
%GHOSTTERMS Gives the matrix elements for the boundary conditions.
% MOET EEN diag, offdiag en bvalue geven voor zowel ux als uy!!!
dom = casedef.dom;

% Checking which boundary the face belongs to
for randID = 1:length(casedef.BC)
    range = dom.getzone(casedef.BC{randID}.zoneID).range;
    if faceIndex >= range(1) && faceIndex <= range(end)
        id =  randID;
        break
    end
end
% Checking which BC applies at that boundary
BC = casedef.BC{id}.velocityKind;
switch BC
    case 'Dirichlet'
        % Determining lambda using the anonymous function
        lambda = getLambda(dom,faceIndex);
        % u
        ughostDiag = 1-lambda; 
        ughostOffdiag = lambda;
        % v
        vghostDiag = 1-lambda; 
        vghostOffdiag = lambda;
        % b values
        uhat = casedef.BC{id}.data.velocity;
        if isa(uhat, 'function_handle')
            pos = dom.fCoord(:,faceIndex);
            b = uhat(pos(1),pos(2));
            ubValue = b(1);
            vbValue = b(2);
        else
            ubValue = uhat(1);
            vbValue = uhat(2);
        end
    case 'Neumann'
        ksi = dom.fXiMag(faceIndex);
        ughostDiag = 1;
        ughostOffdiag = -1;
        vghostDiag = 1;
        vghostOffdiag = -1;
        du = casedef.BC{id}.data.velocity;
        if isa(du, 'function_handle')
            pos = dom.fCoord(:,faceIndex);
            b = du(pos(1),pos(2));
            ubValue = ksi*b(1);
            vbValue = ksi*b(2);
        else
            ubValue = ksi*du(1);
            vbValue = ksi*du(2);
        end
    otherwise
        disp('BC not found');
end
end

