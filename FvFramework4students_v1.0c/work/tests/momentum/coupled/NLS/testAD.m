%% Use automatic differentiation when only first order derivative is needed
func = @(x) NavierStokes(casedef, x);

xAD = myAD(x);
outAD = func(xAD);


% %% Output
% % Function values
% disp('Function values:');
% f
% fAD=getvalue(outAD)
% dfAD=getderivs(outAD)
