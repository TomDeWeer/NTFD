N = 101;
Re = 100;
disp(char("Re="+Reynolds))
disp(char("N="+N))
Nx = N;
Ny = N;
L = 1;
casedef.name = char("Re_"+Reynolds+"_grid_"+Nx+"x"+Ny);
mu = 4;
rho = 10;
Uxtop = Reynolds*mu/(L*rho);
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
plot(yline,uline, 'r-x')
hold on
% paper solution
y = [1.0,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000];
uRe100 = [1.000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717, 0];
plot(y,uRe100,'b-x')
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
    if y==0.5 && x<L && x>0 
        vline = [vline, result.U.data(2,i)/Uxtop];
        xline = [xline, x];
    end
end
figure()
plot(xline,vline,'r-x')
hold on
% paper solution
x = [1.0,0.9688,0.9609,0.9531,0.9453,0.9063,0.8594,0.8047,0.5000,0.2344,0.2266,0.1563,0.0938,0.0781,0.0703,0.0625,0.0000];
vRe100 = [0.00000,-0.05906,-0.07391,-0.08864,-0.10313,-0.16914,-0.22445,-0.24533,0.05454,0.17527,0.17507,0.16077,0.12317,0.10890,0.10091,0.09233,0.00000];
plot(x,vRe100,'b-x')
h = legend('Nonlinear solver','Literature');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
