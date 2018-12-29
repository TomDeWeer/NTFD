function [res,J] = NavierStokes(casedef, x)
% Computes the residuals of the discretized solution to the Navier Stokes
% equations, and the Jacobian.

dom = casedef.dom;

% unpack x
[p, u, v] = getPUV(casedef, x);
Np = length(p);
Nu = length(u);
Nv = length(v);
RC=1;

%% Continuity equations: 
CE = zeros(dom.nPc,1); % contains residuals of equations (1 eq for every cell)
JCE = sparse(double(dom.nPc),double(Np+Nu+Nv)); % J(i,j) is the partial derivative of continuity 
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
    uf = lambda*(ni'*U1) + (1-lambda)*ni'*U2;    
    CE(firstNbC) = CE(firstNbC) + Af*uf;
    JCE(firstNbC, Np+firstNbC) = JCE(firstNbC, Np+firstNbC) + Af*lambda*ni(1); % derivative of u component of U1 to CE of cell 1
    JCE(firstNbC, 2*Np+firstNbC) = JCE(firstNbC, 2*Np+firstNbC)+Af*lambda*ni(2); % derivative of v component of U1 to CE of cell 1
    JCE(firstNbC, Np+secondNbC) = JCE(firstNbC, Np+secondNbC)+ Af*(1-lambda)*ni(1); % derivative of u component of U2 to CE of cell 1
    JCE(firstNbC, 2*Np+secondNbC) = JCE(firstNbC, 2*Np+secondNbC)+ Af*(1-lambda)*ni(2); % derivative of v component of U2 to CE of cell 1
    if secondNbC<=dom.nPc
        % same, but everything is negative because ni is in wrong direction
        CE(secondNbC) = CE(secondNbC) - Af*lambda*ni'*U1 - Af*(1-lambda)*ni'*U2;
        JCE(secondNbC, Np+firstNbC) = JCE(secondNbC, Np+firstNbC) - Af*lambda*ni(1); % derivative of u component of U1 to CE of cell 2
        JCE(secondNbC, 2*Np+firstNbC) = JCE(secondNbC, 2*Np+firstNbC) -Af*lambda*ni(2); % derivative of v component of U1 to CE of cell 2
        JCE(secondNbC, Np+secondNbC) = JCE(secondNbC, Np+secondNbC) -Af*(1-lambda)*ni(1); % derivative of u component of U2 to CE of cell 2
        JCE(secondNbC, 2*Np+secondNbC) = JCE(secondNbC, 2*Np+secondNbC) -Af*(1-lambda)*ni(2); % derivative of v component of U2 to CE of cell 2
    end
    
    if RC && i<dom.nIf
        P1 = p(firstNbC);
        P2 = p(secondNbC);
        % getting further neighbour cells for Rhie Chow
        % getting zeroth neighbour cell
        nb1f = dom.cNbF(4*firstNbC-3:4*firstNbC); 
        nb1fCoords = dom.fCoord(:,nb1f);
        ranking = nb1fCoords'*ni;
        [~, f0] = min(ranking);
        f0 = nb1f(f0);
        Lxi01 = dom.fXiMag(f0);
        possible0Cells = dom.fNbC(2*f0-1:2*f0);
        assert(ismember(firstNbC,possible0Cells))
        if possible0Cells(1) == firstNbC
            zerothNbC = possible0Cells(2);
        else
            zerothNbC = possible0Cells(1);
        end
        % getting third neighbour cell
        nb2f = dom.cNbF(4*secondNbC-3:4*secondNbC); 
        nb2fCoords = dom.fCoord(:,nb2f);
        ranking = nb2fCoords'*ni;
        [~, f3] = max(ranking);
        f3 = nb2f(f3);
        Lxi23 = dom.fXiMag(f3);
        possible3Cells = dom.fNbC(2*f3-1:2*f3);
        assert(ismember(secondNbC,possible3Cells))
        if possible3Cells(1) == secondNbC
            thirdNbC = possible3Cells(2);
        else
            thirdNbC = possible3Cells(1);
        end
        % Applying Rhie Chow
        P0 = p(zerothNbC);
        P1 = p(firstNbC);
        P2 = p(secondNbC);
        P3 = p(thirdNbC);
        Lxi12 = dom.fXiMag(i);
        dpdn = (P2-P1)/(Lxi12);
        dpdn1 = (P2-P0)/(Lxi01 + Lxi12);
        dpdn2 = (P3-P1)/(Lxi12 + Lxi23);
        dpdn_hat = lambda*dpdn1 + (1-lambda)*dpdn2;
        factor = Af;%*omega*(1/af);
        % On first cell
        CE(firstNbC) = CE(firstNbC) + factor*(dpdn -dpdn_hat);
        % Derivative contributions of dpdn_hat
        JCE(firstNbC, zerothNbC) = JCE(firstNbC, zerothNbC) + factor*lambda/(Lxi01 + Lxi12);
        JCE(firstNbC, firstNbC) = JCE(firstNbC, firstNbC) + factor*(1-lambda)/(Lxi12 + Lxi23); 
        JCE(firstNbC, secondNbC) = JCE(firstNbC, secondNbC) - factor*lambda/(Lxi01 + Lxi12); 
        JCE(firstNbC, thirdNbC) = JCE(firstNbC, thirdNbC) - factor*(1-lambda)/(Lxi12 + Lxi23); 
        % Derivative contribution of dpdn
        JCE(firstNbC, firstNbC) = JCE(firstNbC, firstNbC) - factor/Lxi12; 
        JCE(firstNbC, secondNbC)= JCE(firstNbC, secondNbC) + factor/Lxi12; 
        % On second cell
        CE(secondNbC) = CE(secondNbC) - factor*(dpdn -dpdn_hat);
        % Derivative contributions of dpdn_hat
        JCE(secondNbC, zerothNbC) = JCE(secondNbC, zerothNbC) - factor*lambda/(Lxi01 + Lxi12);
        JCE(secondNbC, firstNbC) = JCE(secondNbC, firstNbC) - factor*(1-lambda)/(Lxi12 + Lxi23); 
        JCE(secondNbC, secondNbC) = JCE(secondNbC, secondNbC) + factor*lambda/(Lxi01 + Lxi12); 
        JCE(secondNbC, thirdNbC) = JCE(secondNbC, thirdNbC) + factor*(1-lambda)/(Lxi12 + Lxi23); 
        % Derivative contribution of dpdn
        JCE(secondNbC, firstNbC) = JCE(secondNbC, firstNbC) + factor/Lxi12; 
        JCE(secondNbC, secondNbC)= JCE(secondNbC, secondNbC) - factor/Lxi12; 
    end
end

%% Navier Stokes equations
CNS = zeros(2*dom.nPc,1); % contains residuals for Navier Stokes equations(2 for every physical cell)
% equations are ordered as follows: CNS = [firstcell_u, firstCell_v, secondCell_u, secondCell_v, ...]
JNS = sparse(double(2*dom.nPc), double(Np+Nu+Nv)); % J(i,j) is the partial derivative of NS 
%equation i with respect to x(j)(x = u, v or p)
rho = casedef.material.rho;
% first term: del dot u u 
% syms U1(1) U1(2) U2(1) U2(2) EQ1 EQ2 unSYM;
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
%     unSYM = ni'*(lambda*[U1(1);U1(2)]+(1-lambda)*[U2(1);U2(2)]);
%     EQ1 = Af*(lambda*U1(1)+(1-lambda)*U2(1))*unSYM;
%     EQ2 = Af*(lambda*U1(2)+(1-lambda)*U2(2))*unSYM;
%     jac = jacobian([EQ1,EQ2],[U1(1), U1(2), U2(1), U2(2)]);
%     derivatives = double(subs(jac,[U1(1),U1(2),U2(1),U2(2)],[U1(1),U1(2),U2(1),U2(2)]));
    % derivative of first equation
    JNS(2*firstNbC-1, Np+firstNbC) = JNS(2*firstNbC-1, Np+firstNbC) + ...
        Af*lambda*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) ...
        + ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) ...
        + Af*lambda*ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)); % derivative to U1(1)
    JNS(2*firstNbC-1, Np+Nu+firstNbC) = JNS(2*firstNbC-1,Np+Nu+firstNbC) + ...
        Af*lambda*ni(2)*(U1(1)*lambda - U2(1)*(lambda - 1)); % derivative to U1(2)
    JNS(2*firstNbC-1, Np+secondNbC) = JNS(2*firstNbC-1, Np+secondNbC) ...
        - Af*(lambda - 1)*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) ...
        + ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) ...
        - Af*ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1))*(lambda - 1); % derivative to U2(1)
    JNS(2*firstNbC-1, Np+Nu+secondNbC) = JNS(2*firstNbC-1, Np+Nu+secondNbC) ...
        -Af*ni(2)*(U1(1)*lambda - U2(1)*(lambda - 1))*(lambda - 1); % derivative to U2(2)
    % derivative of second equation
    JNS(2*firstNbC, Np+firstNbC) = JNS(2*firstNbC, Np+firstNbC) ...
        + Af*lambda*ni(1)*(U1(2)*lambda - U2(2)*(lambda - 1)); % derivative to U1(1)
    JNS(2*firstNbC, Np+Nu+firstNbC) = JNS(2*firstNbC, Np+Nu+firstNbC) ...
        + Af*lambda*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) ...
        + ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) ...
        + Af*lambda*ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1)); % derivative to U1(2)
    JNS(2*firstNbC, Np+secondNbC) = JNS(2*firstNbC, Np+secondNbC) ...
        -Af*ni(1)*(U1(2)*lambda - U2(2)*(lambda - 1))*(lambda - 1); % derivative to U2(1)
    JNS(2*firstNbC, Np+Nu+secondNbC) = JNS(2*firstNbC, Np+Nu+secondNbC) ...
        - Af*(lambda - 1)*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) ...
        + ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) ...
        - Af*ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))*(lambda - 1); % derivative to U2(2) 
    if secondNbC<=dom.nPc
        CNS(2*secondNbC-1) = CNS(2*secondNbC-1) - Af*(lambda*U1(1)+(1-lambda)*U2(1))*un;
        CNS(2*secondNbC) = CNS(2*secondNbC) - Af*(lambda*U1(2)+(1-lambda)*U2(2))*un;
%         unSYM = ni'*(lambda*[U1(1);U1(2)]+(1-lambda)*[U2(1);U2(2)]);
%         EQ1 = - Af*(lambda*U1(1)+(1-lambda)*U2(1))*unSYM;
%         EQ2 = - Af*(lambda*U1(2)+(1-lambda)*U2(2))*unSYM;
%         jac = jacobian([EQ1,EQ2],[U1(1), U1(2), U2(1), U2(2)]);
%         derivatives = double(subs(jac,[U1(1),U1(2),U2(1),U2(2)],[U1(1),U1(2),U2(1),U2(2)]));
        % derivative of first equation
        JNS(2*secondNbC-1, Np+firstNbC) = JNS(2*secondNbC-1, Np+firstNbC) ...
            - Af*lambda*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) + ...
            ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) ...
            - Af*lambda*ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)); % derivative to U1(1)
        JNS(2*secondNbC-1, Np+Nu+firstNbC) = JNS(2*secondNbC-1, Np+Nu+firstNbC) ...
            -Af*lambda*ni(2)*(U1(1)*lambda - U2(1)*(lambda - 1)); % derivative to U1(2)
        JNS(2*secondNbC-1, Np+secondNbC) = JNS(2*secondNbC-1, Np+secondNbC) ...
            + Af*(lambda - 1)*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) + ...
            ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) + ...
            Af*ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1))*(lambda - 1); % derivative to U2(1)
        JNS(2*secondNbC-1, Np+Nu+secondNbC) = JNS(2*secondNbC-1, Np+Nu+secondNbC) ...
            + Af*ni(2)*(U1(1)*lambda - U2(1)*(lambda - 1))*(lambda - 1); % derivative to U2(2)
        % derivative of second equation
        JNS(2*secondNbC, Np+firstNbC) = JNS(2*secondNbC, Np+firstNbC) ...
            -Af*lambda*ni(1)*(U1(2)*lambda - U2(2)*(lambda - 1)); % derivative to U1(1)
        JNS(2*secondNbC, Np+Nu+firstNbC) = JNS(2*secondNbC, Np+Nu+firstNbC) ...
            - Af*lambda*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) + ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) - Af*lambda*ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1)); % derivative to U1(2)
        JNS(2*secondNbC, Np+secondNbC) = JNS(2*secondNbC, Np+secondNbC) ...
            + Af*ni(1)*(U1(2)*lambda - U2(2)*(lambda - 1))*(lambda - 1); % derivative to U2(1)
        JNS(2*secondNbC, Np+Nu+secondNbC) = JNS(2*secondNbC, Np+Nu+secondNbC) ...
            + Af*(lambda - 1)*(ni(1)*(U1(1)*lambda - U2(1)*(lambda - 1)) ...
            + ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))) ...
            + Af*ni(2)*(U1(2)*lambda - U2(2)*(lambda - 1))*(lambda - 1); % derivative to U2(2) 
    end
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
    JNS(2*firstNbC-1,firstNbC) = JNS(2*firstNbC-1,firstNbC)  + Af*lambda*ni(1)/rho;
    JNS(2*firstNbC-1,secondNbC) = JNS(2*firstNbC-1,secondNbC)+ Af*(1-lambda)*ni(1)/rho;
    JNS(2*firstNbC,firstNbC) = JNS(2*firstNbC,firstNbC)  + Af*lambda*ni(2)/rho;
    JNS(2*firstNbC,secondNbC) = JNS(2*firstNbC,secondNbC)+ Af*(1-lambda)*ni(2)/rho;
    if secondNbC<=dom.nPc
        CNS(2*secondNbC-1) = CNS(2*secondNbC-1) - Af*(lambda*p1+(1-lambda)*p2)*ni(1)/rho;
        CNS(2*secondNbC) = CNS(2*secondNbC) - Af*(lambda*p1+(1-lambda)*p2)*ni(2)/rho;
        JNS(2*secondNbC-1,firstNbC) = JNS(2*secondNbC-1,firstNbC) -Af*lambda*ni(1)/rho;
        JNS(2*secondNbC-1,secondNbC) = JNS(2*secondNbC-1,secondNbC) -Af*(1-lambda)*ni(1)/rho;
        JNS(2*secondNbC,firstNbC) = JNS(2*secondNbC,firstNbC) -Af*lambda*ni(2)/rho;
        JNS(2*secondNbC,secondNbC) = JNS(2*secondNbC,secondNbC) -Af*(1-lambda)*ni(2)/rho;
    end
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
    % derivative of first equation
    JNS(2*firstNbC-1,Np+firstNbC) = JNS(2*firstNbC-1,Np+firstNbC)  + nu*Af/Lxi; % derivative to U1(1)
    JNS(2*firstNbC-1,Np+secondNbC) = JNS(2*firstNbC-1,Np+secondNbC)  - nu*Af/Lxi; % derivative to U2(1)
    % derivative of second equation
    JNS(2*firstNbC,Np+Nu+firstNbC) = JNS(2*firstNbC,Np+Nu+firstNbC)  + nu*Af/Lxi; % derivative to U1(2)
    JNS(2*firstNbC,Np+Nu+secondNbC) = JNS(2*firstNbC,Np+Nu+secondNbC)  - nu*Af/Lxi; % derivative to U2(2)
    if secondNbC<=dom.nPc
        CNS(2*secondNbC-1) = CNS(2*secondNbC-1) + nu*Af*U2(1)/Lxi - nu*Af*U1(1)/Lxi;
        CNS(2*secondNbC) = CNS(2*secondNbC) + nu*Af*U2(2)/Lxi - nu*Af*U1(2)/Lxi;
        % derivative of first equation
        JNS(2*secondNbC-1,Np+firstNbC) = JNS(2*secondNbC-1,Np+firstNbC)  - nu*Af/Lxi; % derivative to U1(1)
        JNS(2*secondNbC-1,Np+secondNbC) = JNS(2*secondNbC-1,Np+secondNbC)  + nu*Af/Lxi; % derivative to U2(1)
        % derivative of second equation
        JNS(2*secondNbC,Np+Nu+firstNbC) = JNS(2*secondNbC,Np+Nu+firstNbC)  - nu*Af/Lxi; % derivative to U1(2)
        JNS(2*secondNbC,Np+Nu+secondNbC) = JNS(2*secondNbC,Np+Nu+secondNbC) + nu*Af/Lxi; % derivative to U2(2)
    end
end

%% Boundary equations
BCE = []; % boundary condition equations
JBCE = sparse([]); % jacobian of BCE
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
    JBCEi = sparse(double(Np+Nu+Nv),1);% one BCE equation gradient(must be transposed for JBCE)
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
            JBCEi(PC) = lambda;
            JBCEi(GC) = (1-lambda);
            JBCE = [JBCE; JBCEi'];
        case 'Neumann'
            ksi = dom.fXiMag(i);
            fixedPressureGradient = casedef.BC{id}.data.pressure;
            if isa(fixedPressureGradient, 'function_handle')
                pos = dom.fCoord(:,faceIndex);
                fixedPressureGradient = fixedPressureGradient(pos(1),pos(2));
            end
            % Check normalization
            JBCEi = sparse(double(Np+Nu+Nv),1);
            if casedef.BC{id}.isNormalized && i == range(1)
                % apply zero dirichlet BC
%                 lambda = getLambda(dom,i);
%                 BCE = [BCE; lambda*pPC + (1-lambda)*pGC];
%                 JBCEi(PC) = lambda;
%                 JBCEi(GC) = (1-lambda);
%                 JBCE = [JBCE; JBCEi'];
%                 lambda = getLambda(dom,i);
%                 BCE = [BCE; pGC];
%                 JBCEi(GC) = 1;
%                 JBCE = [JBCE; JBCEi'];
                BCE = [BCE; fixedPressureGradient+pPC/(ksi)];
                JBCEi(PC) =  1/ksi;
                JBCE = [JBCE; JBCEi']; 
            else
                BCE = [BCE; fixedPressureGradient-(pGC-pPC)/(ksi)];
                JBCEi(GC) = -1/ksi;
                JBCEi(PC) =  1/ksi;
                JBCE = [JBCE; JBCEi']; 
            end

        otherwise
            disp('BC not found');
    end
    
    BCuv = casedef.BC{id}.velocityKind;
    switch BCuv
        case 'Dirichlet' % fixed 
            % Determining lambda using the anonymous function
            lambda = getLambda(dom,i);
            fixedVelocity = casedef.BC{id}.data.velocity;
            if isa(fixedVelocity, 'function_handle')
                pos = dom.fCoord(:,i);
                fixedVelocity = fixedVelocity(pos(1),pos(2));
            end
            fixedU = fixedVelocity(1); fixedV = fixedVelocity(2);
            BCE = [BCE; lambda*uPC + (1-lambda)*uGC - fixedU];
            JBCEi = sparse(Np+Nu+Nv,1);% one BCE equation gradient(must be transposed for JBCE)
            JBCEi(Np+PC) = lambda;
            JBCEi(Np+GC) = (1-lambda);
            JBCE = [JBCE; JBCEi'];
            BCE = [BCE; lambda*vPC + (1-lambda)*vGC - fixedV];
            JBCEi = sparse(Np+Nu+Nv,1);% one BCE equation gradient(must be transposed for JBCE)
            JBCEi(Np+Nu+PC) = lambda;
            JBCEi(Np+Nu+GC) = (1-lambda);
            JBCE = [JBCE; JBCEi'];
        case 'Neumann'
            ksi = dom.fXiMag(i);
            fixedVelocityGradient = casedef.BC{id}.data.velocity;
            if isa(fixedVelocityGradient, 'function_handle')
                pos = dom.fCoord(:,i);
                fixedVelocityGradient = fixedVelocityGradient(pos(1),pos(2));
            end
            fixedUGradient = fixedVelocityGradient(1);
            fixedVGradient = fixedVelocityGradient(2);
            BCE = [BCE; fixedUGradient-(uGC-uPC)/(ksi)];
            JBCEi = sparse(double(Np+Nu+Nv),1);% one BCE equation gradient(must be transposed for JBCE)
            JBCEi(Np+PC) = 1/ksi;
            JBCEi(Np+GC) = -1/ksi;
            JBCE = [JBCE; JBCEi'];
            BCE = [BCE; fixedVGradient-(vGC-vPC)/(ksi)];
            JBCEi = sparse(Np+Nu+Nv,1);% one BCE equation gradient(must be transposed for JBCE)
            JBCEi(Np+Nu+PC) = 1/ksi;
            JBCEi(Np+Nu+GC) = -1/ksi;
            JBCE = [JBCE; JBCEi'];
        otherwise
            disp('BC not found');
    end
end
res = [CE; CNS; BCE];
J = [JCE; JNS; JBCE];
end
