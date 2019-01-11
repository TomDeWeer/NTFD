%==========================================================================
%
%
%
%==========================================================================
function result = coupledNLS(casedef)
dom = casedef.dom;

% initialize the optimization vector x
% x contains velocities and pressures at every point in the grid
% x = [p ; u; v] where p contains pressures u contains x velocities and y
% contains y velocities

if isfield(casedef,'P') && isfield(casedef,'U')
    uvec = casedef.U.data;
    u = uvec(1,:)';
    v = uvec(2,:)';
    p = casedef.P.data';
    x0 = [p; u; v];
else
    x0 = zeros(3*dom.nC,1);
end
%outputFunc = @(x,optimVals,state) plotFlow(x,casedef, optimVals);
options = optimoptions('fsolve','Display','iter-detailed',...
    'SpecifyObjectiveGradient',true, ...
    'FiniteDifferenceStepSize', 1.e-5, 'TolFun',casedef.iteration.FuncTol, 'OptimalityTolerance', casedef.iteration.OptTol, 'Algorithm','trust-region-dogleg','MaxIterations',5000);%"CheckGradients",true); % 'OutputFcn',outputFunc,
handle = @(x) NavierStokes(casedef, x);
tic
[sol, fval, exitflag, output] = myOwn_fsolve(handle,x0,options);
output.time = toc;
plotFlowFinal(sol,casedef,output);
[p, u, v] = getPUV(casedef, sol);
U = Field(dom.allCells,1);
set(U,[u'; v']);
P = Field(dom.allCells,0);
set(P,p');
result.U=U; result.P=P;
result.sol = sol;
result.output = output;
end
