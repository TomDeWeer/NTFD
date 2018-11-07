function [ ghostDiag, ghostOffdiag, bValue ] ...
    = tGhostTerms( casedef, faceIndex )
%GHOSTTERMS Gives the matrix elements for the boundary conditions.

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
BC = casedef.BC{id}.kind;
switch BC
    case 'Dirichlet'
        % Determining lambda using the anonymous function
        lambda = getLambda(dom,faceIndex);
        ghostDiag = 1-lambda; % IK HEB HET HIER OMGEDRAAID
        ghostOffdiag = lambda;
        bValue = casedef.BC{id}.data.bcval;
    case 'Neumann'
        K = casedef.material.k;
        Af = dom.fArea(faceIndex);
        ksi = dom.fXiMag(faceIndex); % niet nodig
        ghostDiag = K/(ksi);
        ghostOffdiag = -K/(ksi);
        phi = casedef.BC{id}.data.bcval;
        if isa(phi, 'function_handle')
            pos = dom.fCoord(:,faceIndex);
            bValue = phi(pos(1),pos(2));
        else
            bValue = phi;
        end
        
    otherwise
        disp('BC not found');
end


end

