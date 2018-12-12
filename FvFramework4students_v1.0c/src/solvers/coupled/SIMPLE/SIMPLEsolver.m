%==========================================================================
%==========================================================================
function result = SIMPLEsolver(casedef)
dom = casedef.dom;

% Create field objects
% U = Field(dom.allCells,1);	% Velocity [m/s] (2D vector); empty field
% reset(U,[0; 0]); 
% set(U,casedef.U0.data);          % Set to given initial guess

% Create an equation object for holding a vector conservation equation
eqnU = ScalarFvEqn2(dom);
eqnV = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
alpha = casedef.relaxation;
while iterate   
    
    niter = niter+1;
    
    %%% MOMENTUMSOLVER BASED ON GUESSED P %%%
    [Unew, uP, vP] = momguesssolver(casedef);
    set(casedef.U,Unew'); % Put algebraic solution in the Field
    
    %%% CORRECTION FOR P %%%
    Pcorr = pressureCorrSolver(casedef, uP, vP);

    
    %%% UPDATE FIELD VALUES U, V, P %%%
    % update pressures
    set(casedef.P, (casedef.P.data + alpha*Pcorr'));
    % update velocities
    [Uupdated] = updateVelocities(casedef,Pcorr,uP,vP);
    set(casedef.U,Uupdated); % Put algebraic solution in the Field
    
    %%% CHECK CONVERGENCE BASED ON U', V' %%%
    p = casedef.P.data';
    u = casedef.U.data(1,:)';
    v = casedef.U.data(2,:)';
    x = [p; u; v];
    residuals = NavierStokes(casedef, x);
    resnorm = norm(residuals);
    fprintf("It %d : residual norm =  %.12f \n",niter, resnorm)
    if resnorm < casedef.iteration.resTol
        Uconverged = true;
        iterate = false;
        disp("Convergence achieved")
    elseif niter > casedef.iteration.maxniter
        Uconverged = false;
        iterate = false;
        disp("Max iterations reached.")
%    elseif checkstoprequest(stopmon)
%       Uconverged = false;
%       iterate = false;
    end
end % iterate

result.endtime = now; % call datestr(now) for displaying this time 
result.Uconverged = Uconverged;
result.niter = niter;
result.resnorm = resnorm;
result.U = Field(dom.allCells,1);
set(result.U,casedef.U.data);


end



