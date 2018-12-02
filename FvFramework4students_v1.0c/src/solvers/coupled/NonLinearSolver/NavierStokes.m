function [J] = NavierStokes(casedef, x)
% Computes the residuals of the discretized solution to the Navier Stokes
% equations, and the Jacobian.

dom = casedef.dom;

% unpack x
[p, u, v] = getPUV(casedef, x);
Np = length(p);
Nu = length(u);
Nv = length(v);


%% Continuity equations: 
CE = zeros(dom.nPc,1); % contains residuals of equations (1 eq for every cell)
JCE = zeros(dom.nPc,Np+Nu+Nv); % J(i,j) is the partial derivative of continuity 
%equation in cell i with respect to x(j)(x = u, v or p)
for i=1:dom.nF
    Af = dom.fArea(i);
    [firstNbC,secondNbC] = getCells(dom,i);
    lambda = getLambda(dom,i); % points from firstNbc to secondNbc
    U1 = [u(firstNbC); v(firstNbC)];
    U2 = [u(secondNbC); v(secondNbC)];
    % equation for first Nbc
    Xi = dom.fXi(:,i); % pointing correctly for firstNbc
    ni = Xi/norm(Xi,2); % normal
    % TODO: make anonymous function to wrap this up nicely
    CE(firstNbC) = CE(firstNbC) + Af*lambda*(ni'*U1) + Af*(1-lambda)*ni'*U2;
%     JCE(firstNbC, Np+firstNbC) = JCE(firstNbC, Np+firstNbC) + Af*lambda*ni(1); % derivative of u component of U1 to CE of cell 1
%     JCE(firstNbC, 2*Np+firstNbC) = JCE(firstNbC, 2*Np+firstNbC)+Af*lambda*ni(2); % derivative of v component of U1 to CE of cell 1
%     JCE(firstNbC, Np+secondNbC) = JCE(firstNbC, Np+secondNbC)+ Af*(1-lambda)*ni(1); % derivative of u component of U2 to CE of cell 1
%     JCE(firstNbC, 2*Np+secondNbC) = JCE(firstNbC, 2*Np+secondNbC)+ Af*(1-lambda)*ni(2); % derivative of v component of U2 to CE of cell 1
    if secondNbC<=dom.nPc
        % same, but everything is negative because ni is in wrong direction
        CE(secondNbC) = CE(secondNbC) - Af*lambda*ni'*U1 - Af*(1-lambda)*ni'*U2;
        JCE(secondNbC, Np+firstNbC) = JCE(secondNbC, Np+firstNbC) - Af*lambda*ni(1); % derivative of u component of U1 to CE of cell 2
        JCE(secondNbC, 2*Np+firstNbC) = JCE(secondNbC, 2*Np+firstNbC) -Af*lambda*ni(2); % derivative of v component of U1 to CE of cell 2
        JCE(secondNbC, Np+secondNbC) = JCE(secondNbC, Np+secondNbC) -Af*(1-lambda)*ni(1); % derivative of u component of U2 to CE of cell 2
        JCE(secondNbC, 2*Np+secondNbC) = JCE(secondNbC, 2*Np+secondNbC) -Af*(1-lambda)*ni(2); % derivative of v component of U2 to CE of cell 2
    end
end

%% Navier Stokes equations
CNS = zeros(2*dom.nPc,1); % contains residuals for Navier Stokes equations(2 for every physical cell)
% equations are ordered as follows: CNS = [firstcell_u, firstCell_v, secondCell_u, secondCell_v, ...]
JNS = zeros(2*dom.nPc, Np+Nu+Nv); % J(i,j) is the partial derivative of NS 
%equation i with respect to x(j)(x = u, v or p)
rho = casedef.material.rho;
% first term: del dot u u 
for i=1:dom.nF
     Af = dom.fArea(i);
    [firstNbC,secondNbC] = getCells(dom,i);
    lambda = getLambda(dom,i); % points from firstNbc to secondNbc
    U1 = [u(firstNbC); v(firstNbC)];
    U2 = [u(secondNbC); v(secondNbC)];
    % equations for first Nbc
    Xi = dom.fXi(:,i); % pointing correctly for firstNbc
    ni = Xi/norm(Xi,2); % normal
    un = ni'*(lambda*U1+(1-lambda)*U2);
    CNS(2*firstNbC-1) = CNS(2*firstNbC-1)+ Af*(lambda*U1(1)+(1-lambda)*U2(1))*un; % contribution of this face to u equation in firstNbC
    CNS(2*firstNbC) = CNS(2*firstNbC) + Af*(lambda*U1(2)+(1-lambda)*U2(2))*un; % contribution of this face to v equation in firstNbC
%     % derivative of first equation
%     JNS(2*firstNbC-1, Np+2*firstNbC-1) = JNS(2*firstNbC-1, Np+2*firstNbC-1) ...
%         + 2*Af*lambda*lambda*ni(1)*U1(1) + Af*lambda*(ni(1)*(1-lambda)*U2(1)+ni(2)*lambda*U1(2) + ni(2)*(1-lambda)*U2(2)); % derivative to U1(1)
%     JNS(2*firstNbC-1, Np+2*firstNbC) = JNS(2*firstNbC-1, Np+2*firstNbC) + ...
%         Af*(lambda*U1(1)+(1-lambda)*U2(1))*ni(2)*lambda; % derivative to U1(2)
%     JNS(2*firstNbC-1, Np+2*secondNbC-1) = JNS(2*firstNbC-1, Np+2*secondNbC-1) ...
%         + 2*Af*(1-lambda)*(1-lambda)*ni(1); % derivative to U2(1)
%     JNS(2*firstNbC-1, Np+2*secondNbc) = JNS(2*firstNbC-1, Np+2*secondNbc) + ...
%         Af*(lambda*U1(1)+(1-lambda)*U2(1))*ni(2)*(1-lambda); % derivative to U2(2)
%     % derivative of second equation
%     JNS(2*firstNbC, Np+2*firstNbC-1) = JNS(2*firstNbC, Np+2*firstNbC-1) + ...
%         Af*(lambda*U1(2)+(1-lambda)*U2(2))*ni(1)*lambda; % derivative to U1(1)
%     JNS(2*firstNbC, Np+2*firstNbC) = JNS(2*firstNbC, Np+2*firstNbC) + ...
%         Af*lambda*ni(2)*lambda*2*U1(2); % derivative to U1(2)
%     JNS(2*firstNbC, Np+2*secondNbC-1) = JNS(2*firstNbC, Np+2*secondNbC-1) + ...
%         Af*(lambda*U1(2)+(1-lambda)*U2(2))*ni(1)*(1-lambda); % derivative to U2(1)
%     JNS(2*firstNbC, Np+2*secondNbC) = JNS(2*firstNbC, Np+2*secondNbC) + ...
%         2*Af*(1-lambda)*ni(2)*(1-lambda); % derivative to U2(2)    
    if secondNbC<=dom.nPc
        CNS(2*secondNbC-1) = CNS(2*secondNbC-1) - Af*(lambda*U1(1)+(1-lambda)*U2(1))*un;
        CNS(2*secondNbC) = CNS(2*secondNbC) - Af*(lambda*U1(2)+(1-lambda)*U2(2))*un;
    end
    % TODO: jacobiaan invullen voor eerste term
end

% second term: (1/rho) * int(n.p, op rand)
for i=1:dom.nF
     Af = dom.fArea(i);
    [firstNbC,secondNbC] = getCells(dom,i);
    lambda = getLambda(dom,i); % points from firstNbc to secondNbc
    p1 = p(firstNbC);
    p2 = p(secondNbC);
    % equations for first Nbc
    Xi = dom.fXi(:,i); % pointing correctly for firstNbc
    ni = Xi/norm(Xi,2); % normal
    CNS(2*firstNbC-1)= CNS(2*firstNbC-1)+ Af*(lambda*p1+(1-lambda)*p2)*ni(1)/rho; % contribution of this face to u equation in firstNbC
    CNS(2*firstNbC) =  CNS(2*firstNbC) +  Af*(lambda*p1+(1-lambda)*p2)*ni(2)/rho; % contribution of this face to v equation in firstNbC
    if secondNbC<=dom.nPc
        CNS(2*secondNbC-1) = CNS(2*secondNbC-1) - Af*(lambda*p1+(1-lambda)*p2)*ni(1)/rho;
        CNS(2*secondNbC) = CNS(2*secondNbC) - Af*(lambda*p1+(1-lambda)*p2)*ni(2)/rho;
    end
    % TODO: jacobiaan voor tweede term
end

% third term: diffusion
nu = casedef.material.nu;
for i=1:dom.nF
     Af = dom.fArea(i);
    [firstNbC,secondNbC] = getCells(dom,i);
    lambda = getLambda(dom,i); % points from firstNbc to secondNbc
    U1 = [u(firstNbC); v(firstNbC)];
    U2 = [u(secondNbC); v(secondNbC)];
    % equations for first Nbc
    Lxi = dom.fXiMag(i);
    CNS(2*firstNbC-1) = CNS(2*firstNbC-1)+ nu*Af*U1(1)/Lxi - nu*Af*U2(1)/Lxi; % contribution of this face to u equation in firstNbC
    CNS(2*firstNbC) = CNS(2*firstNbC) +    nu*Af*U1(2)/Lxi - nu*Af*U2(2)/Lxi; % contribution of this face to v equation in firstNbC
    if secondNbC<=dom.nPc
        CNS(2*secondNbC-1) = CNS(2*secondNbC-1) + nu*Af*U2(1)/Lxi - nu*Af*U1(1)/Lxi;
        CNS(2*secondNbC) = CNS(2*secondNbC) + nu*Af*U2(2)/Lxi - nu*Af*U1(2)/Lxi;
    end
    % TODO: jacobiaan invullen voor derde term
end

%% Boundary equations
BCE = []; % boundary condition equations
JBCE = []; % jacobian of BCE
% looping over external faces
for i=dom.nIf+1:dom.nF
    [PC,GC] = getCells(dom,i);
    pPC = p(PC);
    pGC = p(GC);
    uPC = u(PC);
    uGC = u(GC);
    vPC = v(PC);
    vGC = v(GC);
    % Checking which boundary the face belongs to
    for randID = 1:length(casedef.BC)
        range = dom.getzone(casedef.BC{randID}.zoneID).range;
        if i >= range(1) && i <= range(end)
            id =  randID;
            break
        end
    end
    % Checking which BC applies at that boundary
    % Pressure
    BCp = casedef.BC{id}.pressureKind;
    JBCEi = zeros(Np+Nu+Nv,1);% one BCE equation gradient(must be transposed for JBCE)
    switch BCp
        case 'Dirichlet'
            % Determining lambda using the anonymous function
            lambda = getLambda(dom,i);
            fixedPressure = casedef.BC{id}.data.pressure;
            if isa(fixedPressure, 'function_handle')
                pos = dom.fCoord(:,i);
                fixedPressure = fixedPressure(pos(1),pos(2));
            end
            BCE = [BCE; lambda*pPC + (1-lambda)*pGC - fixedPressure];
        case 'Neumann'
            ksi = dom.fXiMag(i);
            fixedPressureGradient = casedef.BC{id}.data.pressure;
            if isa(fixedPressureGradient, 'function_handle')
                pos = dom.fCoord(:,faceIndex);
                fixedPressureGradient = fixedPressureGradient(pos(1),pos(2));
            end
            BCE = [BCE; fixedPressureGradient-(pGC-pPC)/(ksi)];
        otherwise
            disp('BC not found');
    end
    BCuv = casedef.BC{id}.velocityKind;
    JBCEi = zeros(Np+Nu+Nv,1);% one BCE equation gradient(must be transposed for JBCE)
    switch BCuv
        case 'Dirichlet' % fixed pressure
            % Determining lambda using the anonymous function
            lambda = getLambda(dom,i);
            fixedVelocity = casedef.BC{id}.data.velocity;
            if isa(fixedVelocity, 'function_handle')
                pos = dom.fCoord(:,i);
                fixedVelocity = fixedVelocity(pos(1),pos(2));
            end
            fixedU = fixedVelocity(1); fixedV = fixedVelocity(2);
            BCE = [BCE; lambda*uPC + (1-lambda)*uGC - fixedU];
            BCE = [BCE; lambda*vPC + (1-lambda)*vGC - fixedV];
        case 'Neumann'
            ksi = dom.fXiMag(i);
            fixedVelocityGradient = casedef.BC{id}.data.velocity;
            if isa(fixedVelocityGradient, 'function_handle')
                pos = dom.fCoord(:,faceIndex);
                fixedVelocityGradient = fixedVelocityGradient(pos(1),pos(2));
            end
            fixedUGradient = fixedVelocityGradient(1);
            fixedVGradient = fixedVelocityGradient(2);
            BCE = [BCE; fixedUGradient-(uGC-uPC)/(ksi)];
            BCE = [BCE; fixedVGradient-(vGC-vPC)/(ksi)];
        otherwise
            disp('BC not found');
    end
    
end
J = [CE; CNS; BCE];
end
