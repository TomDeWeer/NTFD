clear all; close all;
path = [pwd '\Figures\Couette'];
[~, ~, ~] = mkdir(path);

% Create a mesh
Lx = 1;
Ly = 1;
Nx = 20;
Ny = 20;
p0 = 1;
Uxtop = 0; % m/s
mu = 1;
rho = 1;
nu = mu/rho;
dx = Lx/Nx;
dy = Ly/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');
casedef.name = char("NLS_Couette_grid_"+Nx+"x"+Ny);

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[Uxtop*ones(1,U.elcountzone);zeros(1,U.elcountzone)]);
casedef.U = U; % initial guess
% Define material properties
casedef.material.nu = nu;  % viscosity [dynamic???]
casedef.material.rho = rho; % density [kg/m^3]

jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
% casedef.BC{jBC}.kind   = 'Dirichlet';
% casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
casedef.BC{jBC}.velocityKind   = 'Neumann';
casedef.BC{jBC}.data.velocity = [0 , 0]; % uniform velocity inlet
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.velocityKind   = 'Neumann';
casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
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
casedef.BC{jBC}.data.velocity = [Uxtop , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;

% Set up iteration parameters
casedef.iteration.regularization = 1.e-5;
casedef.iteration.restol     = 1.e-5;

% Call solver
result = coupledNLS(casedef);


% Plot result
% ux
figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.2;
fvmplotfield(result.U,scale,lw, 1);
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
colorbar('TickLabelInterpreter', 'latex');
% saveas(gcf,fullfile(path,'Couette_50x50.png'));


% 
% 
dPdx = -p0/Lx;
realU = @(y) (-dPdx/(2*mu))*y*(Ly-y)+Uxtop*y/Ly; % Analytic solution
realP = @(x) p0 + x*dPdx;


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

% convergence experiment results
Ns = [2, 4, 8, 16, 32];
MaximumErrors = [];
AverageErrors = [];
for N = Ns
    fprintf('N: %f \n',N)
    Lx = 1;
    Ly = 1;
    Nx = N;
    Ny = N;
    
    dx = Lx/Nx;
    dy = Ly/Ny;
    seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
    casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
    mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
    % Create domain from mesh
    casedef.dom = newdomain(mesh,'MyDomain');
    casedef.name = char("NLS_Couette_grid_"+Nx+"x"+Ny);

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


    jBC = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'WESTRAND';
    % casedef.BC{jBC}.kind   = 'Dirichlet';
    % casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = [0 , 0];
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = p0; 
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'OOSTRAND';
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
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
    casedef.BC{jBC}.data.velocity = [Uxtop , 0]; % no slip condition
    casedef.BC{jBC}.pressureKind   = 'Neumann';
    casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
    casedef.BC{jBC}.isNormalized = 0;

    % Set up iteration parameters
    casedef.iteration.regularization = 1.e-5;
    casedef.iteration.restol     = 1.e-12;


    % Call solver
    result = coupledNLS(casedef);


    maxErr = 0;
    avgErr = 0;
    Ni=0;
    PErr = Field(casedef.dom.allCells, 0);
    set(PErr,zeros(1,PErr.elcountzone));
    comp_perr = zeros(1, size(PErr.data,2));
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
            comp_perr(i) = abs(result.P.data(i) - realP(x));
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
    set(PErr,comp_perr);
    figure; hold on; axis image; colormap(jet(50));
    scale = 'lin'; lw = 0.2; colorbar(); title('Pressure Error')
    fvmplotfield(PErr,scale,lw);
    
    
    figure; hold on; axis image; colormap(jet(50));
    scale = 'lin'; lw = 0.2;
    fvmplotfield(result.U,scale,lw, 2);
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    colorbar('TickLabelInterpreter', 'latex');
    
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
% saveas(gcf,fullfile(path,'Couette_error_convergence.png'));