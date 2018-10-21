function [ ghostDiag, ghostOffdiag, bValue ] ...
    = ghostTerms( casedef, faceIndex )
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
        ghostDiag = lambda;
        ghostOffdiag = 1-lambda;
    case 'Neumann'
        ghostDiag = -1/dom.fXiMag(faceIndex);
        ghostOffdiag = 1/dom.fXiMag(faceIndex);
    otherwise
        disp('BC not found');
end
bValue = casedef.BC{id}.data.bcval;

end

