close all;
N = 51;
Re = 0.01;
disp(char("Re="+Re))
disp(char("N="+N))
Nx = N;
Ny = N;
L = 1;
casedef.name = char("Re_"+Re+"_grid_"+Nx+"x"+Ny);
mu = 4;
rho = 10;
Uxtop = Re*mu/(L*rho);
% Reynolds = Uxtop*L*rho/mu;
%fprintf("Reynolds number: %.2f \n",Reynolds)
nu = mu/rho;
dx = L/Nx;
dy = L/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 L],Ny,0.99,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[Uxtop*ones(1,U.elcountzone);zeros(1,U.elcountzone)]);
casedef.U = U; % initial guess

% Define material properties
casedef.material.nu = nu;  % viscosity [dynamic???]
casedef.material.rho = rho; % density [kg/m^3]



% Define boundary conditions
% There are boundary conditions for velocity and for pressure
% Wall: dirichlet velocity and neumann pressure in n
% Prescribed inlet: dirichlet velocity and dirichlet pressure
% Unprescribed outlet: neumann velocity in n and dirichlet pressure
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
% casedef.BC{jBC}.kind   = 'Dirichlet';
% casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = @(x,y) [0 , 0]; % uniform velocity inlet
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0 , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 1;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [Uxtop , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
% Set up iteration parameters
casedef.iteration.FuncTol     = 1.e-8;
casedef.iteration.OptTol      = 1.e-6;
casedef.iteration.regularization = 1.e-5;
result = coupledNLS(casedef);

% compare with data from paper 
% u velocity at x=1/2
uline = [];
yline = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if x==0.5 && y<L && y>0 
        uline = [uline, result.U.data(1,i)/Uxtop];
        yline = [yline, y];
    end
end
figure()
plot(yline,uline, 'b')
hold on
% paper solution
ylit = [0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375];
uRe001 = [-3.85e-2,-6.96e-2,-9.69e-2,-1.23e-1,-1.47e-1,-1.71e-1,-1.92e-1,-2.05e-1,-2.06e-1,-1.86e-1,-1.32e-1,-3.24e-2,1.27e-1,3.55e-1,6.51e-1];
plot(ylit,uRe001,'k:x')
% v velocity at y=1/2
h = legend('Nonlinear solver','Literature','Interpreter','Latex');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');


% v velocity at y=1/2
vline = [];
xline = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if y< 0.5+dy/2 && y>0.5-dy/2 && x<L && x>0 
        vline = [vline, result.U.data(2,i)/Uxtop];
        xline = [xline, x];
    end
end
figure()
plot(xline,vline,'b')
hold on
% paper solution
xlit = [0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375];
vRe001 = [9.46e-2,1.56e-1,1.83e-1,1.79e-1,1.52e-1,1.09e-1,5.66e-2,6.37e-6,-5.66e-2,-1.09e-1,-1.52e-1,-1.79e-1,-1.83e-1,-1.56e-1,-9.46e-2];
plot(xlit,vRe001,'k:x')
h = legend('Nonlinear solver','Literature');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');