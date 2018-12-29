function [step,quadObj,normStep,normStepScal, normdNewton] = ...
    myOwn_dogleg(nvar,F,JAC,grad,Delta,scalMat,reg)
%DOGLEG approximately solves trust region subproblem via a dogleg approach.
%
%   DOGLEG finds an approximate solution d to the problem:
%
%     min_d      f + g'd + 0.5*d'Bd
%
%     subject to ||Dd|| <= Delta
%
%   where g is the gradient of f, B is a Hessian approximation of f, D is a
%   diagonal scaling matrix and Delta is a given trust region radius.

%   Copyright 1990-2012 The MathWorks, Inc.

% NOTE: The scaling matrix D above is called scalMat in this routine.

% Compute scaled gradient and other scaled terms.
gradscal = grad./scalMat;
gradscal2 = gradscal./scalMat;
normgradscal = norm(gradscal);
normdNewton = 0;
if normgradscal >= eps
    % First compute the Cauchy step (in scaled space).
    dCauchy = -(Delta/normgradscal)*gradscal;
    JACvec = JAC*gradscal2;
    denom = Delta*(JACvec'*JACvec);
    tauterm = normgradscal^3/denom;
    tauC = min(1,tauterm);
    dCauchy = tauC*dCauchy;
    
    % Compute quadratic objective at Cauchy point.
    JACvec = JAC*(dCauchy./scalMat); 
    objCauchy = gradscal'*dCauchy + 0.5*(JACvec'*JACvec);   
    normdCauchy = min(norm(dCauchy),Delta);
else
    % Set Cauchy step to zero step and continue.
    objCauchy = 0;
    normdCauchy = 0;
    dCauchy = zeros(nvar,1);
end

if Delta - normdCauchy < eps;
    % Take the Cauchy step if it is at the boundary of the trust region.
    step = dCauchy; quadObj = objCauchy;
else
    % Compute the Gauss-Newton step (in scaled space).
    % DEES WAS DUS ECHT SLECHTE CODE
    % Disable the warnings about conditioning for singular and
    % nearly singular matrices
%     warningstate1 = warning('off','MATLAB:nearlySingularMatrix');
%     warningstate2 = warning('off','MATLAB:singularMatrix');
%     warningstate3 = warning('off','MATLAB:rankDeficientMatrix');
    % dNewton = -JAC\F; 
    A = JAC'*JAC + sparse(reg*speye(size(JAC)));
    b = -JAC'*F;
%     A = JAC;
%     b = -F;
    % preconditioning
%     setup.type = 'nofill';
%     [L,U] = ilu(JAC'*JAC,setup);
%     %Minv = inv(U)*inv(L);
%     Aprec = sparse(size(A,1), size(A,2));
%     for colindex=1:size(A,2)
%         col = A(:,colindex);
%         Aprec(:,colindex) = U\(sparse(L\col));
%     end
%     bprec = U\(L\b);
%     dNewton = (Aprec)\(bprec);
%     tol = 1e-5;
%     maxit = 1000;
%     restart = [];
%     setup.type = 'ilutp';
%     setup.udiag = 1;
%     setup.droptol = 1.e-2;
    %[L,U] = ilu(A,setup);
    % dNewton = gmres(A,b,restart,tol,maxit, A');
    if condest(JAC)<condest(A)
        dNewton = JAC\-F;
    else
        dNewton = A\b;
    end
    
     % voila se, vele beter
    % Restore the warning states to their original settings
    normdNewton = norm(dNewton);
%     fprintf("Estimated condition of Jacobian: %.3e \n",condest(JAC))
%     fprintf("Estimated condition of LS-Jacobian: %.3e \n",condest(A))
%     fprintf("Norm of dNewton: %.3f \n",normdNewton)
%     fprintf("gmres relative residual norm: %.3f \n",norm(A*dNewton -b)/norm(b))

%     fprintf("Real relative residual norm: %.3f \n",norm(JAC*dNewton + F)/norm(F))
%     fprintf("Regularization parameter: %.3e \n",reg)

    
%     warning(warningstate1)
%     warning(warningstate2)
%     warning(warningstate3)
    
    dNewton = dNewton.*scalMat;     % scale the step
    
    if any(~isfinite(dNewton))
        % Take the Cauchy step if the Gauss-Newton step gives bad values.
        step = dCauchy; quadObj = objCauchy;
    else 
        normdNewt = norm(dNewton);
        if normdNewt <= Delta
            % Use the Newton direction as the trial step
            step = dNewton;
        else
            % Find the intersect point along dogleg path.
            Delta2 = Delta^2;
            normdCauchy2 = min(normdCauchy^2,Delta2);
            normdNewt2 = normdNewt^2;
            dCdN = dCauchy'*dNewton;
            dCdNdist2 = max((normdCauchy2+normdNewt2-2*dCdN),0);
            
            if dCdNdist2 == 0
                tauI = 0;
            else
                % Stable method for solving 1-D quadratic
                a = 0.5*dCdNdist2;
                b = dCdN - normdCauchy2;
                c = 0.5*(normdCauchy2 - Delta2);
                q = -0.5*(b + sign(b)*sqrt(b^2 - 4*a*c));
                if b > 0
                    tauI = c/q; 
                else
                    tauI = q/a;
                end

                % Make sure we take a finite step,
                % (check for poorly scaled/infinite directions).
                if ~isfinite(tauI)
                    tauI = 0; % Take Cauchy step
                end
            end
            step = dCauchy + tauI*(dNewton-dCauchy);
        end
        % Compute quadratic objective at trial point.
        JACvec = JAC*(step./scalMat);
        quadObj = gradscal'*step + 0.5*(JACvec'*JACvec);
        
        % Compare Cauchy step and trial step (Newton or intersection)
        if objCauchy < quadObj
            step = dCauchy; quadObj = objCauchy;
        end
    end
end

% The step computed was the scaled step.  Unscale it.
normStepScal = norm(step);
step = step./scalMat;
normStep = norm(step);




