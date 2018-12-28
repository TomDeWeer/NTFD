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
% Analyzes the influence of the convergence criterium on the performance of
% the SIMPLE algorithm. The errors of the solutions and the computing time 
% needed are plotted for various NS-based tolerances.
% All analyses are performed using the following parameters:
%     - 10x10 and 20x20 grid
%     - Lx, Ly = 1m
%     - mu = 4Ns/m2
%     - rho = 10kg/m
%     - dpdx = -0.5Pa/m
%     - p1 = 5Pa
%     - no scaling, orthogonal grids
%     - alpha = 0.1
%     - dt = 0.01
%     - Max. iterations = 1000 (unlimited)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General commands
clc;
close all;
clear;
% mkdir('..\Figuren\Opgave 2.2');

%% convergence experiment results
Ns = [16, 32, 64];
Tolpowers = [1, 2, 3, 4, 5, 6, 7, 8];
MaximumErrors = zeros(length(Ns),length(Tolpowers));
AverageErrors = zeros(length(Ns),length(Tolpowers));
TotalTime = zeros(length(Ns),length(Tolpowers));
for i = 1:length(Ns)
    for j = 1:length(Tolpowers)
        N = Ns(i);
        tol = Tolpowers(j);
        fprintf('N: \n',N,'\n Tol: \n 10^-',tol)
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
        for k=1:casedef.dom.nC
            coord = casedef.dom.cCoord(:,k);
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
        casedef.iteration.maxniter = 1000;
        casedef.iteration.resTol = 10^(-tol);
        casedef.iteration.dt = 0.01;
        % relaxation factor
        casedef.relaxation = 0.1;

        % Call solver
        tic
        result = SIMPLEsolver(casedef);
        time = toc;
        TotalTime(i,j) = time;
        fprintf("Time: %.3f \n",time)

        maxErr = 0;
        avgErr = 0;
        Ni=0;
        Err = Field(casedef.dom.allCells, 0);
        set(Err,zeros(1,U.elcountzone));
        comp_err = zeros(1, size(Err.data,2));
        for k=1:result.U.dom.nC
            x = result.U.dom.cCoord(1,k);
            y = result.U.dom.cCoord(2,k);
            % Only keep the interior cells
            if x>0 && y>0 && x<Lx && y<Ly
                Uxapprox = result.U.data(1,k);
                Uxexact = realU(y);
                % Compute relative error
                err = abs((Uxexact-Uxapprox));
                comp_err(k) = err/abs(Uxexact);
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

        avgErr = avgErr/Ni;
        fprintf("Maximum error: %.10f \n", maxErr)
        MaximumErrors(i,j) = maxErr;
        AverageErrors(i,j) = avgErr;
    end
end


figure()
hold on
loglog(-0.1,-0.1,'k','LineWidth',1,'displayname','16 by 16 grid');
loglog(-0.1,-0.1,'b','LineWidth',1,'displayname','32 by 32 grid');
loglog(-0.1,-0.1,'color',[0 128/255 0],'LineWidth',1,'displayname','64 by 64 grid');
loglog(10.^(-Tolpowers),AverageErrors(1,:),'Marker','.','LineStyle','-','color','k','HandleVisibility','off')
loglog(10.^(-Tolpowers),AverageErrors(2,:),'Marker','.','LineStyle','-','color','b','HandleVisibility','off')
loglog(10.^(-Tolpowers),AverageErrors(3,:),'Marker','.','LineStyle','-','color',[0 128/255 0],'HandleVisibility','off')
xlabel('Tolerance on convergence criterium [-]','Interpreter','latex');
ylabel('Error [m/s]','Interpreter','latex');
legend({'16 by 16 grid','32 by 32 grid','64 by 64 grid'},'Interpreter','latex','Location','eastoutside');
title('Error','interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XDir','reverse','XScale', 'log','YScale', 'log')
% saveas(gcf,fullfile(path,'Couette_error_convergence.png'));

figure()
hold on
loglog(-0.1,-0.1,'k','LineWidth',1,'displayname','16 by 16 grid');
loglog(-0.1,-0.1,'b','LineWidth',1,'displayname','32 by 32 grid');
loglog(-0.1,-0.1,'color',[0 128/255 0],'LineWidth',1,'displayname','64 by 64 grid');
loglog(10.^(-Tolpowers),MaximumErrors(1,:),'-x','color','k','HandleVisibility','off')
loglog(10.^(-Tolpowers),MaximumErrors(2,:),'-x','color','b','HandleVisibility','off')
loglog(10.^(-Tolpowers),MaximumErrors(3,:),'-x','color',[0 128/255 0],'HandleVisibility','off')
xlabel('Tolerance on convergence criterium [-]','Interpreter','latex');
ylabel('Error [m/s]','Interpreter','latex');
legend({'16 by 16 grid','32 by 32 grid','64 by 64 grid'},'Interpreter','latex','Location','eastoutside');
title('Error','interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XDir','reverse','XScale', 'log','YScale', 'log')
% saveas(gcf,fullfile(path,'Couette_error_convergence.png'));

figure()
hold on
loglog(10.^(-Tolpowers),TotalTime(1,:),'-k.')
loglog(10.^(-Tolpowers),TotalTime(2,:),'-b.')
xlabel('Tolerance on convergence criterium [-]','Interpreter','latex');
ylabel('Time [s]','Interpreter','latex');
legend('10 by 10 grid','20 by 20 grid', 'Interpreter','latex')
title('Total computing time','interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XDir','reverse','XScale', 'log','YScale', 'log')
% saveas(gcf,fullfile(path,'Couette_error_convergence.png'));
