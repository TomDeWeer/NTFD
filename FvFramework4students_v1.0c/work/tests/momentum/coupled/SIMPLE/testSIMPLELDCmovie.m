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
% Visualises the convergence of the SIMPLE algorithm for LDC flow and
% compares it to results from literature.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

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
Uxtop = Re*mu/(L*rho);
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
set(U,[Uxtop*ones(1,U.elcountzone);zeros(1,U.elcountzone)]);
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
casedef.BC{jBC}.isNormalized = true;
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
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = @(x,y) [Uxtop , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;

% compare with data from paper 
% paper solution
% u velocity at x=1/2
y = [1.0,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000];
uRe100 = [1.000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717, 0];
% v velocity at y=1/2
x = [1.0,0.9688,0.9609,0.9531,0.9453,0.9063,0.8594,0.8047,0.5000,0.2344,0.2266,0.1563,0.0938,0.0781,0.0703,0.0625,0.0000];
vRe100 = [0.00000,-0.05906,-0.07391,-0.08864,-0.10313,-0.16914,-0.22445,-0.24533,0.05454,0.17527,0.17507,0.16077,0.12317,0.10890,0.10091,0.09233,0.00000];

figure()
subplot(2,2,1);     % Subplot 1
hold on
plot(y,uRe100,'k:x')
plot(0,0,':k');
title('$U_x$ at $x=0.5$','interpreter','latex')
subplot(2,2,2);     % Subplot 2
hold on
plot(x,vRe100,'k:x')
plot([0 1],[0 0],'b');
title('$U_y$ at $y=0.5$','interpreter','latex')
subplot(2,2,3);     % Subplot 3
hold on
plot(0,0,':k');
title('$P$ at $x=0.5$','interpreter','latex')
subplot(2,2,4);     % Subplot 4
hold on
plot([0 1],[0 0],'b');
title('$P$ at $y=0.5$','interpreter','latex')
for iter=1:400
    
    % Set up iteration parameters
    casedef.iteration.maxniter = 0;
    casedef.iteration.resTol = 1.e-6;
    casedef.iteration.dt = 0.01;
    % relaxation factor
    casedef.relaxation = 0.1;
    % Call solver
    result = SIMPLEsolver(casedef);
    fprintf("It %d : residual norm =  %.12f \n",iter, result.resnorm)
    
    uline = [];
    pline = [];
    yline = [];
    for i=1:result.U.dom.nC
        x = result.U.dom.cCoord(1,i);
        y = result.U.dom.cCoord(2,i);
        % Only keep the interior cells
        if x>0.497 && x < 0.503 && y<L && y>0 
            uline = [uline, result.U.data(1,i)/Uxtop];
            pline = [pline, result.P.data(i)];
            yline = [yline, y];
        end
    end
    subplot(2,2,1);     % Subplot 1
    children = get(gca, 'children');
    delete(children(1));
    plot(yline,uline, 'b')
    h = legend('Literature','SIMPLE');
    set(h,'interpreter','Latex','FontSize',11)
    set(gca,'TickLabelInterpreter', 'latex');
    subplot(2,2,3);     % Subplot 3
    children = get(gca, 'children');
    delete(children(1));
    plot(yline,pline, 'b')
    set(gca,'TickLabelInterpreter', 'latex');

    vline = [];
    pline = [];
    xline = [];
    for i=1:result.U.dom.nC
        x = result.U.dom.cCoord(1,i);
        y = result.U.dom.cCoord(2,i);
        % Only keep the interior cells
        if y>0.499 && y < 0.501 && x<L && x>0 
            vline = [vline, result.U.data(2,i)/Uxtop];
            pline = [pline, result.P.data(i)];
            xline = [xline, x];
        end
    end
    subplot(2,2,2);     % Subplot 2
    children = get(gca, 'children');
    delete(children(1));
    plot(xline,vline,'b')
    h = legend('Literature','SIMPLE');
    set(h,'interpreter','Latex','FontSize',11)
    set(gca,'TickLabelInterpreter', 'latex');
    subplot(2,2,4);     % Subplot 4
    children = get(gca, 'children');
    delete(children(1));
    plot(xline,pline,'b')
    set(gca,'TickLabelInterpreter', 'latex');
    
    pause(0.05)
end
