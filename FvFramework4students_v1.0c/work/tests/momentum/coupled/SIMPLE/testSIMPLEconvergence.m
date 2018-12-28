%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      NUMERICAL TECHNIQUES IN FLUID DYNAMICS       %%%%%%%%%%%%
%%%%%%%%%%%%             Koen Devesse, Tom De Weer             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SIMPLE algorithm                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test case: Grid convergence
%
% Analyzes the grid convergence of the SIMPLE algorithm. The errors of the
% solutions and the computing time needed are plotted.
% All analyses are performed using the following parameters:
%     - Nx, Ny varying between 2 and 64.
%     - Lx, Ly = 1m
%     - mu = 4Ns/m2
%     - rho = 10kg/m
%     - dpdx = -0.5Pa/m
%     - p1 = 5Pa
%     - no scaling, orthogonal grids
%     - alpha = 0.1
%     - dt = 0.01
%     - Convergence criterium = 10^-6
%     - Max. iterations = 500
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General commands
clc;
close all;
clear;
% mkdir('..\Figuren\Opgave 2.2');

%% convergence experiment results
Ns = [2, 4, 8, 16, 32];
MaximumErrors = [];
MaximumErrorsPos = [];
AverageErrors = [];
TotalTime = [];
StepTime = [];
for N = Ns
    fprintf('N: %f \n',N)
    Lx = 1;
    Ly = 1;
    Nx = N;
    Ny = N;
    Uxtop = 0; % m/s
    mu = 4;
    p1 = 5;
    p2 = 0;
    dPdx = (p2-p1)/Lx;
    rho = 10;
    nu = mu/rho;
    dx = Lx/Nx;
    dy = Ly/Ny;
    realU = @(y) (-dPdx/(2*mu))*y*(Ly-y); % Analytic solution
    seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
    casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
    mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
    % Create domain from mesh
    casedef.dom = newdomain(mesh,'MyDomain');

    % Set up initial fields
    U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
    set(U,[zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
    casedef.U = U; % initial guess
    P = Field(casedef.dom.allCells,0); % Pressure
    P0 = [];
    for i=1:casedef.dom.nC
        coord = casedef.dom.cCoord(:,i);
        x = coord(1);
        y = coord(2);
    %     P0 = [P0, pfunc(x,y)];
    %     P0 = [P0, pfunc(x,y) + sin(pi*x/dx )];
    %     P0 = [P0, sin(pi*x/dx )];
        P0 = [P0, 0];
    end
    set(P,P0)
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
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = @(x,y) [0 , 0]; % uniform velocity inlet
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = p1; % prescribed pressure
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'OOSTRAND';
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = p2; % prescribed pressure
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
    casedef.BC{jBC}.data.velocity = [0 , 0]; % no slip condition
    casedef.BC{jBC}.pressureKind   = 'Neumann';
    casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
    casedef.BC{jBC}.isNormalized = 0;

    % Set up iteration parameters
    casedef.iteration.maxniter = 500;
    casedef.iteration.resTol = 1.e-6;
    casedef.iteration.dt = 0.01;
    % relaxation factor
    casedef.relaxation = 0.1;

    % Call solver
    tic
    result = SIMPLEsolver(casedef);
    time = toc;
    TotalTime = [TotalTime, time];
    StepTime = [StepTime, time/result.niter];
    fprintf("Time: %.3f \n",time)


    maxErr = 0;
    avgErr = 0;
    Ni=0;
    Err = Field(casedef.dom.allCells, 0);
    set(Err,zeros(1,U.elcountzone));
    comp_err = zeros(1, size(Err.data,2));
    for i=1:result.U.dom.nC
        x = result.U.dom.cCoord(1,i);
        y = result.U.dom.cCoord(2,i);
        % Only keep the interior cells
        if x>0 && y>0 && x<Lx && y<Ly
            Uxapprox = result.U.data(1,i);
            Uxexact = realU(y);
            % Compute relative error
            err = abs((Uxexact-Uxapprox));
            comp_err(i) = err/abs(Uxexact);
            % Compute average error
            avgErr = avgErr + err;
            Ni = Ni+1;
            % Update maximum error
            if err>maxErr
                maxErr = err;
                maxErrx = x;
                maxErry = y;
            end
        end
    end
    set(Err,comp_err);
    figure; hold on; axis image; colormap(jet(50));
    scale = 'lin'; lw = 0.2; colorbar(); title('Error')
    fvmplotfield(Err,scale,lw);
    
    avgErr = avgErr/Ni;
    fprintf("Maximum error: %.10f \n", maxErr)
    MaximumErrors = [MaximumErrors, maxErr];
    MaximumErrorsPos = [MaximumErrorsPos, [maxErrx; maxErry]];
    AverageErrors = [AverageErrors, avgErr];
end


figure()
loglog(Ns,MaximumErrors,':x','color','b')
hold on
loglog(Ns,AverageErrors,'-b.')
xlabel('N [-]','Interpreter','latex');
ylabel('Error [m/s]','Interpreter','latex');
legend('Maximum','Average', 'Interpreter','latex')
title('Error','interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
% saveas(gcf,fullfile(path,'Couette_error_convergence.png'));

figure()
hold on
loglog(Ns,TotalTime,'-b.')
loglog(Ns,StepTime,':x','color','b')
xlabel('N [-]','Interpreter','latex');
ylabel('Time [s]','Interpreter','latex');
legend('Total','Per iteration', 'Interpreter','latex')
title('Computing time','interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XScale', 'log','YScale', 'log')
% saveas(gcf,fullfile(path,'Couette_error_convergence.png'));
