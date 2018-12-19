clear all; close all;
path = [pwd '\Figuren\Couette'];
[~, ~, ~] = mkdir(path);
% Create a mesh
Lx = 1;
Ly = 1;
Nx = 20;
Ny = 20;
dPdx = -10;
Uxtop = 2; % m/s
mu = 4;
rho = 10;
nu = mu/rho;
dx = Lx/Nx;
dy = Ly/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[ones(1,U.elcountzone);ones(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros
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


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
% casedef.BC{jBC}.kind   = 'Dirichlet';
% casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = [0, 0];
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = [0 , 0];
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = [0, 0];
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = [Uxtop, 0];

% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.UTol     = 1.e-12;
casedef.iteration.dt = 100;

% Call solver
result = momentumsolver(casedef);


% Plot result
% ux
figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.2;
fvmplotfield(result.U,scale,lw, 1);
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
colorbar('TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Couette_50x50.png'));



% figure()
% quiver(casedef.dom.cCoord(1,:),...
%     casedef.dom.cCoord(2,:),result.U.data(1,:),result.U.data(2,:));

% comparing to analytical solution
% line = [];
% lineyloc = [];
% % Loop over all cells
% for i=1:result.U.dom.nC
%     x = result.U.dom.cCoord(1,i);
%     y = result.U.dom.cCoord(2,i);
%     % Select the cells at the vertical line
%     if x < (0.5*Lx+dx) && x>(0.5*Lx) && y>0 && y<Ly
%         line = [line, result.U.data(1,i)];
%         % Store cell x coordinate
%         lineyloc = [lineyloc , y];
%     end
% end
% % Plot results
% figure()
% hold on;
% plot(lineyloc, line)
% xlabel("y")
% ylabel("u_x")
% title("Fully developed flow profile")
% 
% 
realU = @(y) (-dPdx/(2*mu))*y*(Ly-y)+Uxtop*y/Ly; % Analytic solution
% realLine = zeros(size(line));
% for i=1:length(line)
%     realLine(i) = realU(lineyloc(i));
% end
% plot(lineyloc, realLine)


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
        err = abs((Uxexact-Uxapprox)/Uxexact);
        comp_err(i) = Uxexact-Uxapprox;
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
avgErr = avgErr/Ni;
fprintf("Maximum error: %.10f \n", maxErr)

figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.2; colorbar(); title('Error')
fvmplotfield(Err,scale,lw);

% 
% % comparing to analytical solution: horizontal line
% line = zeros(Nx,1);
% linexloc = zeros(Nx,1);
% % Loop over all cells
% for i=1:result.U.dom.nC
%     x = result.U.dom.cCoord(1,i);
%     y = result.U.dom.cCoord(2,i);
%     % Select the cells at the vertical line
%     if y < (0.5*Ly+dy) && y>(0.5*Ly) && x>0 && x<Lx
%         ix = round((x+0.5*dx)/dx);
%         % Store temperature
%         line(ix) = result.U.data(1,i);
%         % Store cell x coordinate
%         linexloc(ix) = x;
%     end
% end
% % Plot results
% figure()
% hold on;
% plot(linexloc, line)
% xlabel("y")
% ylabel("u_x")
% title("Horizontal line")
% 
% 
% 


% convergence experiment results
Ns = [2, 4, 8, 16, 32, 64];
MaximumErrors = [];
AverageErrors = [];
for N = Ns
    fprintf('N: %f \n',N)
    Lx = 1;
    Ly = 1;
    Nx = N;
    Ny = N;
    dPdx = -10;
    Uxtop = 2; % m/s
    mu = 4;
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
    % reset(U,[1;0.2]); 
    % reset(U,[0; 0]);                          % Reset with all zeros
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


    % Define boundary conditions
    jBC = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'WESTRAND';
    % casedef.BC{jBC}.kind   = 'Dirichlet';
    % casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
    casedef.BC{jBC}.kind   = 'Neumann';
    casedef.BC{jBC}.data.bcval = [0, 0];
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'OOSTRAND';
    casedef.BC{jBC}.kind   = 'Neumann';
    casedef.BC{jBC}.data.bcval = [0 , 0];
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'ZUIDRAND';
    casedef.BC{jBC}.kind   = 'Dirichlet';
    casedef.BC{jBC}.data.bcval = [0, 0];
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'NOORDRAND';
    casedef.BC{jBC}.kind   = 'Dirichlet';
    casedef.BC{jBC}.data.bcval = [Uxtop, 0];

    % Set up iteration parameters
    casedef.iteration.maxniter = 1000;
    casedef.iteration.UTol     = 1.e-8;
    casedef.iteration.dt = 20;

    % Call solver
    result = momentumsolver(casedef);


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
            comp_err(i) = err;
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
    AverageErrors = [AverageErrors, avgErr];
end


figure()
loglog(Ns,MaximumErrors, "kx")
hold on
loglog(Ns,AverageErrors, "rx")
xlabel('N [-]','Interpreter','latex');
ylabel('Error [K]','Interpreter','latex');
legend('Maximum','Average', 'Interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Couette_error_convergence.png'));

