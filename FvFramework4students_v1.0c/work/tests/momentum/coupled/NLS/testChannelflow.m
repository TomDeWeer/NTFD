clear all; close all;
maxUxErrors = [];
maxUyErrors = [];
maxPErrors = [];
Ni = [5, 10, 20, 40, 80, 160];
times = [];
for N=Ni
    % Create a mesh
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

    % Call solver
    % result = momentumsolver(casedef);

    % Plot result
    % ux
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Ux"); colorbar();
%     fvmplotfield(result.U,scale,lw, 1);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Uy"); colorbar();
%     fvmplotfield(result.U,scale,lw, 2);


    result = coupledNLS(casedef);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 0.2; colorbar();
%     fvmplotfield(result.P,scale,lw);
%     set(gca,'TickLabelInterpreter', 'latex');
%     xlabel('x [m]','Interpreter','latex');
%     ylabel('y [m]','Interpreter','latex');
%     colorbar('TickLabelInterpreter', 'latex');
    
    times = [times, result.output.time];
    
    res = NavierStokes(casedef,result.sol);
    fprintf("Final residual value ((res'*res)): %.10f \n", (res'*res))
    % X-Velocity error
    UXErr = Field(casedef.dom.allCells, 0);
    set(UXErr,zeros(1,UXErr.elcountzone));
    computed_Uxerr = zeros(1, size(UXErr.data,2));
    for i=1:result.U.dom.nC
        x = result.U.dom.cCoord(1,i);
        y = result.U.dom.cCoord(2,i);
        % Only keep the interior cells
        if x>0 && y>0 && x<Lx && y<Ly
            Uxapprox = result.U.data(1,i);
            Uxexact = ufunc(x,y);
            % Compute relative error
            computed_Uxerr(i) = abs(Uxexact-Uxapprox);
        end
    end
    set(UXErr,computed_Uxerr);
    avgUXErr = mean2(computed_Uxerr);
    maxUXErr = max(max(computed_Uxerr));
    fprintf("Maximum Ux error: %.10f \n", maxUXErr)
    fprintf("Average Ux error: %.10f \n", avgUXErr)
    maxUxErrors = [maxUxErrors, maxUXErr];

    % Y-Velocity error
    UYErr = Field(casedef.dom.allCells, 0);
    set(UYErr,zeros(1,UYErr.elcountzone));
    computed_Uyerr = zeros(1, size(UYErr.data,2));
    for i=1:result.U.dom.nC
        x = result.U.dom.cCoord(1,i);
        y = result.U.dom.cCoord(2,i);
        % Only keep the interior cells
        if x>0 && y>0 && x<Lx && y<Ly
            Uyapprox = result.U.data(2,i);
            Uyexact = vfunc(x,y);
            % Compute relative error
            computed_Uyerr(i) = abs(Uxexact-Uxapprox);
        end
    end
    set(UYErr,computed_Uyerr);
    avgUYErr = mean2(computed_Uyerr);
    maxUYErr = max(max(computed_Uyerr));
    fprintf("Maximum Uy error: %.10f \n", maxUYErr)
    fprintf("Average Uy error: %.10f \n", avgUYErr)
    maxUyErrors = [maxUyErrors, maxUYErr];
    % Pressure error
    PErr = Field(casedef.dom.allCells, 0);
    set(PErr,zeros(1,PErr.elcountzone));
    computed_Perr = zeros(1, size(PErr.data,2));
    for i=1:result.U.dom.nC
        x = result.U.dom.cCoord(1,i);
        y = result.U.dom.cCoord(2,i);
        % Only keep the interior cells
        if x>0 && y>0 && x<Lx && y<Ly
            Papprox = result.P.data(i);
            Pexact = pfunc(x,y);
            % Compute relative error
            computed_Perr(i) = abs(Pexact-Papprox);
        end
    end
    set(PErr,computed_Perr);
    avgPErr = mean2(computed_Perr);
    maxPErr = max(max(computed_Perr));
    fprintf("Maximum P error: %.10f \n", maxPErr)
    fprintf("Average P error: %.10f \n", avgPErr)

    maxPErrors = [maxPErrors, maxPErr];
    figure; hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("Ux error"); colorbar();
    fvmplotfield(UXErr,scale,lw);
    figure; hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("Uy error"); colorbar();
    fvmplotfield(UYErr,scale,lw);
    figure; hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("P error"); colorbar();
    fvmplotfield(PErr,scale,lw);
%     xi = [];
%     pi = [];
%     for i=1:length(result.P.data)
%         p = result.P.data(i);
%         coord = casedef.dom.cCoord(:,i);
%         x = coord(1);
%         y = coord(2);
%         if y < 0.5 + dx && y>0.5
%             if x>0
%                 xi = [xi, x];
%                 pi = [pi, p];
%             else%     xi = [];
%     pi = [];
%     for i=1:length(result.P.data)
%         p = result.P.data(i);
%         coord = casedef.dom.cCoord(:,i);
%         x = coord(1);
%         y = coord(2);
%         if y < 0.5 + dx && y>0.5
%             if x>0
%                 xi = [xi, x];
%                 pi = [pi, p];
%             else
%                 xi = [x, xi];
%                 pi = [p, pi];
%             end
%         end
%     end
%     figure;
%     plot(xi, pi)

%                 xi = [x, xi];
%                 pi = [p, pi];
%             end
%         end
%     end
%     figure;
%     plot(xi, pi)
end

figure()
loglog(Ni, maxUxErrors, 'bx')
hold on
loglog(Ni,maxUyErrors, 'rx')
loglog(Ni,maxPErrors, 'kx')
xlabel('N','Interpreter','latex' )
h = legend('max($\mid U_x - \bar{U}_x \mid$)','max($\mid U_y - \bar{U}_y \mid$)','max($\mid P - \bar{P} \mid$)');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
figure()
plot(Ni,times,'k-x')
xlabel('N','Interpreter', 'Latex')
ylabel('Computation time [s]', 'Interpreter', 'Latex')
set(gca,'TickLabelInterpreter', 'latex');

