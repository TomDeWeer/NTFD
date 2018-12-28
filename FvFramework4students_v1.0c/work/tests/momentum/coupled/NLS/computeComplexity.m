% compare timings for LDC and NLS
% Create a mesh
Nis = [2, 4, 8, ];
Res = [ 1, 2, 4];
for Reynolds = Res
    iterations = [];
    times = [];
    for N = Nis
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
        casedef.iteration.OptTol      = 1.e-8;
        casedef.iteration.regularization = 1.e-5;
        result = coupledNLS(casedef);
        times = [times, result.output.time];
        iterations = [iterations, result.output.funcCount];
    end
    figure()
    plot(Nis, times)
    figure()
    plot(Nis, iterations)
    figure()
    plot(Nis, times./iterations)
end
    
    


