clear all; close all;

% Create a mesh
Lx = 1;
Ly = 1;
Nx = 5;
Ny = 5;
mu = 4;
dPdx = -10;
Uxtop = 0;
p0 = 5;
rho = 10;
nu = mu/rho;
dx = Lx/Nx;
dy = Ly/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[ones(1,U.elcountzone);ones(1,U.elcountzone)]);
casedef.U = U; % initial guess
gradP = Field(casedef.dom.allCells,1); % Pressure
gradP0 = [];
for i=1:casedef.dom.nC
    gradP0 = [gradP0, [dPdx; 0]];
end
set(gradP,gradP0)
casedef.gradP = gradP;
% Define material properties
casedef.material.nu = nu;  % viscosity [dynamic???]
casedef.material.rho = rho; % density [kg/m^3]

% actual solution
pfunc = @(x,y) p0 + dPdx*x;
ufunc = @(x,y) (-dPdx/(2*mu))*y*(Ly-y)+Uxtop*y/Ly;
vfunc = @(x,y) 0;

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
casedef.BC{jBC}.data.velocity = @(x,y) [ufunc(x,y) , vfunc(x,y)]; % uniform velocity inlet
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = pfunc; % prescribed pressure
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.velocityKind   = 'Neumann';
casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = pfunc; % prescribed pressure
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0 , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0 , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative

% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.UTol     = 1.e-10;
casedef.iteration.dt = 20;

% Call solver
% result = momentumsolver(casedef);

% Plot result
% ux
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("Ux"); colorbar();
% fvmplotfield(result.U,scale,lw, 1);
% figure; hold on; axis off; axis equal; colormap(jet(50));
% scale = 'lin'; lw = 1; title("Uy"); colorbar();
% fvmplotfield(result.U,scale,lw, 2);




p = zeros(casedef.dom.nC,1);
u = zeros(casedef.dom.nC,1);
v = zeros(casedef.dom.nC,1);
for i=1:casedef.dom.nC
    x = casedef.dom.cCoord(1,i);
    y = casedef.dom.cCoord(2,i);
    p(i) = pfunc(x,y);
    u(i) = ufunc(x,y);
    v(i) = vfunc(x,y);
end

x = [p; u; v];
[CE, CNS, BCE] = NavierStokes(casedef, x);
fprintf("CE norm: %.10f \n",norm(CE))
fprintf("NS norm: %.10f \n", norm(CNS))
fprintf("BCE norm: %.10f \n", norm(BCE))
