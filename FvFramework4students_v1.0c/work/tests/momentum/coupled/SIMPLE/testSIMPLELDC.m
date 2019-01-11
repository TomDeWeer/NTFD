%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      NUMERICAL TECHNIQUES IN FLUID DYNAMICS       %%%%%%%%%%%%
%%%%%%%%%%%%             Koen Devesse, Tom De Weer             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SIMPLE algorithm                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test case: Lid Driven Cavity
%
% Uses the SIMPLE algorithm for the LDC problem.
% Two cases are studied: Re=0.01(51x51grid) and Re=100(51x51grid). The
% results are compared to values from literature.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General commands
clc;
close all;
clear;
% mkdir('..\Figuren\Opgave 2.2');

%% Re = 0.01
N = 51;
Re = 0.01;
disp(char("Re="+Re))
disp(char("N="+N))
Nx = N;
Ny = N;
L = 1;
casedef.name = char("Re_"+Re+"_grid_"+Nx+"x"+Ny);
mu = 100;
rho = 1;
Uxtop = 1;
% Reynolds = Uxtop*L*rho/mu;
%fprintf("Reynolds number: %.2f \n",Reynolds)
nu = mu/rho;
dx = L/Nx;
dy = L/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 L],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[Uxtop*zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
casedef.U = U; % initial guess
P = Field(casedef.dom.allCells,0); % Pressure
set(P,zeros(1,P.elcountzone));
casedef.P = P;

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
casedef.BC{jBC}.data.velocity = @(x,y) [Uxtop , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
% Set up iteration parameters
casedef.iteration.maxniter = 800;
casedef.iteration.resTol = 1.e-7;
casedef.iteration.dt = 1e-5;
% relaxation factor
casedef.relaxation = 0.01;

% Call solver
tic
result = SIMPLEsolver(casedef);
fprintf("Time: %.3f \n",toc)




%%% compare with data from paper %%%
% paper solution %
% u velocity at x=1/2
ylit = [0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375];
uRe001 = [-3.85e-2,-6.96e-2,-9.69e-2,-1.23e-1,-1.47e-1,-1.71e-1,-1.92e-1,-2.05e-1,-2.06e-1,-1.86e-1,-1.32e-1,-3.24e-2,1.27e-1,3.55e-1,6.51e-1];
% v velocity at y=1/2
xlit = [0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375];
vRe001 = [9.46e-2,1.56e-1,1.83e-1,1.79e-1,1.52e-1,1.09e-1,5.66e-2,6.37e-6,-5.66e-2,-1.09e-1,-1.52e-1,-1.79e-1,-1.83e-1,-1.56e-1,-9.46e-2];
% SIMPLE solution %
uline = [];
yline = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if x>0.499 && x < 0.501 && y<L && y>0 
        uline = [uline, result.U.data(1,i)/Uxtop];
        yline = [yline, y];
    end
end
vline = [];
xline = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if y>0.499 && y < 0.501 && x<L && x>0 
        vline = [vline, result.U.data(2,i)/Uxtop];
        xline = [xline, x];
    end
end
% Figure %
figure()
subplot(1,2,1);     % Subplot 1
hold on
plot(ylit,uRe001,'k:x')
plot(yline,uline, 'b')
h = legend('Literature','SIMPLE');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
title('$U_x$ at $x=0.5$','interpreter','latex')
subplot(1,2,2);     % Subplot 2
hold on
plot(xlit,vRe001,'k:x')
plot(xline,vline,'b')
h = legend('Literature','SIMPLE');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
title('$U_y$ at $y=0.5$','interpreter','latex')
% Plot field
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Ux"); colorbar();
fvmplotfield(result.U,scale,lw, 1);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Uy"); colorbar();
fvmplotfield(result.U,scale,lw, 2);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("P"); colorbar();
fvmplotfield(result.P,scale,lw);

%% Re = 100
N = 51;
Re = 100;
disp(char("Re="+Re))
disp(char("N="+N))
Nx = N;
Ny = N;
L = 1;
casedef.name = char("Re_"+Re+"_grid_"+Nx+"x"+Ny);
mu = 0.01;
rho = 1;
Uxtop = 1;
% Reynolds = Uxtop*L*rho/mu;
%fprintf("Reynolds number: %.2f \n",Reynolds)
nu = mu/rho;
dx = L/Nx;
dy = L/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 L],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[Uxtop*zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
casedef.U = U; % initial guess
P = Field(casedef.dom.allCells,0); % Pressure
set(P,zeros(1,P.elcountzone));
casedef.P = P;

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
casedef.BC{jBC}.data.velocity = @(x,y) [Uxtop , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
% Set up iteration parameters
casedef.iteration.maxniter = 800;
casedef.iteration.resTol = 1.e-10;
casedef.iteration.dt = 0.05;
% relaxation factor
casedef.relaxation = 0.1;

% Call solver
tic
result = SIMPLEsolver(casedef);
fprintf("Time: %.3f \n",toc)


%%% compare with data from paper %%%
% paper solution %
% u velocity at x=1/2
ylit = [1.0,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000];
uRe100 = [1.000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717, 0];
% v velocity at y=1/2
xlit = [1.0,0.9688,0.9609,0.9531,0.9453,0.9063,0.8594,0.8047,0.5000,0.2344,0.2266,0.1563,0.0938,0.0781,0.0703,0.0625,0.0000];
vRe100 = [0.00000,-0.05906,-0.07391,-0.08864,-0.10313,-0.16914,-0.22445,-0.24533,0.05454,0.17527,0.17507,0.16077,0.12317,0.10890,0.10091,0.09233,0.00000];
% SIMPLE solution %
uline = [];
yline = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if x>0.499 && x < 0.501 && y<L && y>0 
        uline = [uline, result.U.data(1,i)/Uxtop];
        yline = [yline, y];
    end
end
vline = [];
xline = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if y>0.499 && y < 0.501 && x<L && x>0 
        vline = [vline, result.U.data(2,i)/Uxtop];
        xline = [xline, x];
    end
end
% Figure %
figure()
subplot(1,2,1);     % Subplot 1
hold on
plot(ylit,uRe100,'k:x')
plot(yline,uline, 'b')
h = legend('Literature','SIMPLE');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
title('$U_x$ at $x=0.5$','interpreter','latex')
subplot(1,2,2);     % Subplot 2

hold on
plot(xlit,vRe100,'k:x')
plot(xline,vline,'b')
h = legend('Literature','SIMPLE');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
title('$U_y$ at $y=0.5$','interpreter','latex')
% Plot field
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Ux"); colorbar();
fvmplotfield(result.U,scale,lw, 1);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Uy"); colorbar();
fvmplotfield(result.U,scale,lw, 2);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("P"); colorbar();
fvmplotfield(result.P,scale,lw);
