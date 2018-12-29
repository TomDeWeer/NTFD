function [uGhost, vGhost] = UcorrGhostTerms(casedef, faceIndex, uPhys, vPhys)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dom = casedef.dom;

% Checking which boundary the face belongs to
for randID = 1:length(casedef.BC)
    range = dom.getzone(casedef.BC{randID}.zoneID).range;
    if faceIndex >= range(1) && faceIndex <= range(end)
        id =  randID;
        break
    end
end
BC = casedef.BC{id}.velocityKind;
switch BC
    case 'Dirichlet'
        % Determining lambda using the anonymous function
        lambda = getLambda(dom,faceIndex);
        % Prescribed values
        uhat = casedef.BC{id}.data.velocity;
        if isa(uhat, 'function_handle')
            pos = dom.fCoord(:,faceIndex);
            b = uhat(pos(1),pos(2));
            uValue = b(1);
            vValue = b(2);
        else
            uValue = uhat(1);
            vValue = uhat(2);
        end
        % u
        uGhost = (uValue - lambda*uPhys) / (1 - lambda);
        % v
        vGhost = (vValue - lambda*vPhys) / (1 - lambda);
    case 'Neumann'
        ksi = dom.fXiMag(faceIndex);
        du = casedef.BC{id}.data.velocity;
        if isa(du, 'function_handle')
            pos = dom.fCoord(:,faceIndex);
            du_described = du(pos(1),pos(2));
        else
            du_described = du;
        end
        % u
        uGhost = -ksi*du_described(1) + uPhys;
        % v
        vGhost = -ksi*du_described(2) + vPhys;
    otherwise
        disp('BC not found');
end
end

