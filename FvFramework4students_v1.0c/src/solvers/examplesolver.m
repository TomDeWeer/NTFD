%==========================================================================
%
% Example solver using the FVMLab framework 4 students
%
% Purpose: Provides code structure for solving a scalar conservation
%          equation using data structures provided by the framework.
%
% by Frederik Rogiers
%
%==========================================================================
function result = examplesolver(casedef)

dom = casedef.dom;

% Create field objects
T = Field(dom.allCells,0);      % Temperature [K] (scalar); empty field
reset(T,0);                     % Reset with all zeros

% ... Create all other required data structures

% Create an equation object for holding a scalar conservation equation
eqn = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
while iterate   
   
   niter = niter+1;
   
   % Set all terms to zero
   reset(eqn); 
   nC = dom.nC;
   nIf = dom.nIf;
   nBf = dom.nBf;
   Adiag = zeros(nC, 1);
   AoffdiagI = zeros(2*nIf, 1);
   AoffdiagB = zeros(2*nBf, 1);
   b = zeros(nC,1);
   kappa = casedef.material.k;
   % Compute coefficients for physical cell eqns and add them to eqn object
   for i= 1:nIf
       % Calculating terms of the equations
       Af = dom.fArea(i);
       Lxi = dom.fXiMag(i);
       anb = kappa*Af/Lxi;
       % Placing terms in matrix
       % Diagonal
       firstNbC = dom.fNbC(2*i-1);
       secondNbC = dom.fNbC(2*i);
       Adiag(firstNbC) = Adiag(firstNbC) - anb;
       Adiag(secondNbC) = Adiag(secondNbC) - anb;
       % Offdiagonal
       AoffdiagI(2*i-1) = anb;
       AoffdiagI(2*i) = anb;
   end
   % Compute coefficients for ghost cell eqns and add them to eqn object
   for i = 1:nBf
       % Calculating terms of the equations
       Af = dom.fArea(nIf+i);
       Lxi = dom.fXiMag(nIf+i);
       anb = kappa*Af/Lxi;
       % Iterating through the neighbouring cells
       firstNbC = dom.fNbC(2*(nIf+i)-1);
       secondNbC = dom.fNbC(2*(nIf+i));
       for indexCell = [firstNbC, secondNbC]
           % Checking whether it's a physical cell or a ghost cell
           if indexCell <= dom.nPc
                Adiag(indexCell) = Adiag(indexCell) - anb;
           else % If it's a ghost cell
               % Checking which boundary the face belongs to
               for randID = 1:length(casedef.BC)
                   range = casedef.dom.getzone(casedef.BC{randID}.zoneID).range;
                   if nIf+i >= range(1) && nIf+i <= range(end)
                       id =  randID;
                       break
                   end
               end
               % Checking which BC applies at that boundary
               BC = casedef.BC{id}.kind;
               switch BC
                   case 'Dirichlet'
                       lambda = 1 - norm(dom.cCoord(:,indexCell) ...
                           - dom.fCoord(:,nIf+i)) / dom.fXiMag(nIf+i);
                       Adiag(indexCell) = lambda;
                       GCvglOffdiag = 1-lambda;
                   case 'Neumann'
                       Adiag(indexCell) = -1/dom.fXiMag(nIf+i);
                       GCvglOffdiag = 1/dom.fXiMag(nIf+i);
                   otherwise
                       disp('BC not found');
               end
               b(indexCell) = casedef.BC{id}.data.bcval;
           end
       end
       % Filling in the offdiagonal elements
       AoffdiagB(2*i-1) = anb;
       AoffdiagB(2*i) = GCvglOffdiag;
   end
   A = [Adiag; AoffdiagI; AoffdiagB];
   
   eqn.adata = A;
   eqn.bdata = b;
   
   % Create a matlab sparse linear system from the eqn object
   [A,b] = to_msparse(eqn);
   x = get(T);
   x = x';
   
   % Check tolerance and iteration count
   TRes = b-A*x;
   TResnorm = norm(TRes);         
   if TResnorm < casedef.iteration.TTol
      Tconverged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      Tconverged = false;
      iterate = false;
%    elseif checkstoprequest(stopmon)
%       Tconverged = false;
%       iterate = false;
   else
      x = A\b; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(T,x'); % Put algebraic solution in the Field
   end
      
   
end % iterate

result.endtime = now; % call datestr(now) for displaying this time 
result.Tconverged = Tconverged;
result.niter = niter;
result.TResnorm = TResnorm;
result.TRes = Field(dom.allCells,0);
set(result.TRes,TRes');
result.T = T;


end



