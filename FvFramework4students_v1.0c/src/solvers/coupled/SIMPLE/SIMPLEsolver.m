%==========================================================================
%==========================================================================
function result = SIMPLEsolver(casedef)
dom = casedef.dom;

% Create field objects
% U = Field(dom.allCells,1);	% Velocity [m/s] (2D vector); empty field
% reset(U,[0; 0]); 
% set(U,casedef.U0.data);          % Set to given initial guess

% Create an equation object for holding a vector conservation equation
eqnU = ScalarFvEqn2(dom);
eqnV = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
alpha = casedef.relaxation;

% figure()
% subplot(1,2,1);     % Subplot 1
% axis([-0.1 1.1 -0.5 5.5]);
% axis square
% hold on
% plot([0 1],[5 0],':k');
% plot(0,0,':k');
% plot(1,1,':k');
% title('Pressure in x','interpreter','latex')
% subplot(1,2,2);     % Subplot 2
% axis([0 1 0 0.2]);
% axis square
% hold on
% plot(0:0.01:1,(5/(2*4))*(0:0.01:1).*(1-(0:0.01:1)),'k');
% plot([0 1],[0 0],'b');
% title('$U_x$','interpreter','latex')

while iterate   
    
    niter = niter+1;
    
    %%% MOMENTUMSOLVER BASED ON GUESSED P %%%
    [Unew, uP, vP] = momguesssolver(casedef);
    set(casedef.U,Unew'); % Put algebraic solution in the Field
    
    %%% CORRECTION FOR P %%%
    Pcorr = pressureCorrSolver(casedef, uP, vP);

    
    %%% UPDATE FIELD VALUES U, V, P %%%
    % update pressures
    set(casedef.P, (casedef.P.data + alpha*Pcorr'));
    % update velocities
    [Uupdated] = updateVelocities(casedef,Pcorr,uP,vP);
%     Ucorrfield = Field(dom.allCells,1);
%     set(Ucorrfield,Uupdated-Unew');
%     PcorrField = Field(dom.allCells,0);
%     set(PcorrField, Pcorr')
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Pcorr"); colorbar();
%     fvmplotfield(PcorrField,scale,lw);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Ux correction"); colorbar();
%     fvmplotfield(Ucorrfield,scale,lw, 1);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Ux"); colorbar();
%     fvmplotfield(casedef.U,scale,lw, 1);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Uy correction"); colorbar();
%     fvmplotfield(Ucorrfield,scale,lw, 2);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Uy"); colorbar();
%     fvmplotfield(casedef.U,scale,lw, 2);
    
    set(casedef.U,Uupdated); % Put algebraic solution in the Field
    
%     dx = casedef.dom.cCoord(2,1);
%     xi = [];
%     pi = [];
%     ValueAtMiddle = [];
%     for i=1:length(casedef.P.data)
%         p = casedef.P.data(i);
%         coord = casedef.dom.cCoord(:,i);
%         x = coord(1);
%         y = coord(2);
%         if y == 0.5 + dx
%             if x>0
%                 xi = [xi, x];
%                 pi = [pi, p];
%             else
%                 xi = [x, xi];
%                 pi = [p, pi];
%             end
%         end
%         if x == 0.5 + dx
%             if y>0
%                 ValueAtMiddle = [ValueAtMiddle, [y; casedef.U.data(1,i)]];
%             else
%                 ValueAtMiddle = [[y; casedef.U.data(1,i)], ValueAtMiddle];
%             end
%         end
%     end
%     subplot(1,2,1);     % Subplot 1
%     children = get(gca, 'children');
%     delete(children(1));
%     plot(xi, pi, 'b')
%     subplot(1,2,2);     % Subplot 2
%     children = get(gca, 'children');
%     delete(children(1));
%     plot(ValueAtMiddle(1,:),ValueAtMiddle(2,:),'b');
%     pause(0.05)
    
%     PcorrField = Field(dom.allCells,0);
%     set(PcorrField, Pcorr')
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("Pcorr"); colorbar();
%     fvmplotfield(PcorrField,scale,lw);
%     figure; hold on; axis off; axis equal; colormap(jet(50));
%     scale = 'lin'; lw = 1; title("P"); colorbar();
%     fvmplotfield(casedef.P,scale,lw);
    
    %%% CHECK CONVERGENCE BASED ON U', V' %%%
    p = casedef.P.data';
    u = casedef.U.data(1,:)';
    v = casedef.U.data(2,:)';
    x = [p; u; v];
    residuals = NavierStokes(casedef, x);
    resnorm = norm(residuals);
%     fprintf("It %d : residual norm =  %.12f \n",niter, resnorm)
    if resnorm < casedef.iteration.resTol
        Uconverged = true;
        iterate = false;
        disp("Convergence achieved")
    elseif niter > casedef.iteration.maxniter
        Uconverged = false;
        iterate = false;
        disp("Max iterations reached.")
%    elseif checkstoprequest(stopmon)
%       Uconverged = false;
%       iterate = false;
    end
end % iterate

result.endtime = now; % call datestr(now) for displaying this time 
result.Uconverged = Uconverged;
result.niter = niter;
result.resnorm = resnorm;
result.U = Field(dom.allCells,1);
set(result.U,casedef.U.data);
result.P = Field(dom.allCells,0);
set(result.P, casedef.P.data)

end



