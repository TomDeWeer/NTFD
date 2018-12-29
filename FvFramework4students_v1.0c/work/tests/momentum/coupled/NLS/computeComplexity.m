% compare timings for LDC and NLS
% Create a mesh
Nis = [2, 4, 8, ];
Res = [ 1, 2, 4];

NLStimes = [];
SIMPLEtimes = [];
for N = Nis
    
    Lx = 1;
    Ly = 1;
    Nx = N;
    Ny = N;
    mu = 4;
    dPdx = -10;
    Uxtop = 2;
    p0 = 5;
    rho = 10;
    nu = mu/rho;
    dx = Lx/Nx;
    dy = Ly/Ny;
    seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
    casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
    casedef.name = char("Channel__grid_"+Nx+"x"+Ny);

    mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
    % Create domain from mesh
    casedef.dom = newdomain(mesh,'MyDomain');

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
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = [0, 0]; % uniform velocity inlet
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = pfunc; % prescribed pressure
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'OOSTRAND';
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = pfunc; % prescribed pressure
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'ZUIDRAND';
    casedef.BC{jBC}.velocityKind   = 'Dirichlet';
    casedef.BC{jBC}.data.velocity = @(x,y) [ufunc(x,y), vfunc(x,y)]; % no slip condition
    casedef.BC{jBC}.pressureKind   = 'Neumann';
    casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'NOORDRAND';
    casedef.BC{jBC}.velocityKind   = 'Dirichlet';
    casedef.BC{jBC}.data.velocity = @(x,y) [ufunc(x,y), vfunc(x,y)]; % no slip condition
    casedef.BC{jBC}.pressureKind   = 'Neumann';
    casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
    casedef.BC{jBC}.isNormalized = 0;
    % Set up iteration parameters
    casedef.iteration.FuncTol     = 1.e-8;
    casedef.iteration.OptTol      = 1.e-6;
    casedef.iteration.regularization = 1.e-5;
    result = coupledNLS(casedef);
    NLStimes = [times, result.output.time];
    
    clear casedef;
    Lx = 1;
    Ly = 1;
    Nx = N;
    Ny = N;
    mu = 4;
    dPdx = -10;
    Uxtop = 2;
    p0 = 5;
    rho = 10;
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
    casedef.iteration.maxniter = 500;
    casedef.iteration.resTol = 1.e-6;
    casedef.iteration.dt = 0.01;
    % relaxation factor
    casedef.relaxation = 0.1;

    % Call solver
    tic
    result = SIMPLEsolver(casedef);
    time = toc;
    
end
figure()
plot(Nis, times)
figure()
plot(Nis, iterations)
figure()
plot(Nis, times./iterations)

    
    


