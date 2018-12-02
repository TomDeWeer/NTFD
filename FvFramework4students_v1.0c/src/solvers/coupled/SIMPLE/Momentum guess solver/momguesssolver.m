function [Unew,uP,vP] = momguesssolver(casedef)
dom = casedef.dom;

% Create field objects
U = Field(dom.allCells,1);	% Velocity [m/s] (2D vector); empty field
reset(U,[0; 0]); 
set(U,casedef.U.data);      % Set to given initial guess

% Create an equation object for holding a scalar conservation equation
eqnU = ScalarFvEqn2(dom);
eqnV = ScalarFvEqn2(dom);


% Set all terms to zero
reset(eqnU); 
reset(eqnV); 

[Au, bu, Av, bv] = mgMatrixMaker(casedef);

% Create a matlab sparse linear system from the eqn object
eqnU.adata = Au;
[Au,~] = to_msparse(eqnU);
eqnV.adata = Av;
[Av,~] = to_msparse(eqnV);

Un = casedef.U.data;
u = Un(1,:)';
v = Un(2,:)';

[Cu,Ru] = qr(Au,bu);
unew = Ru\Cu;
[Cv,Rv] = qr(Av,bv);
vnew = Rv\Cv;
Unew = [unew, vnew];
uP = diag(Au);
vP = diag(Av);
end
