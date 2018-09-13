% augmented Lagrangian method to minimise f(x) under the linear constraint Ax=b
% f needs to be a matlab function that takes a Matlab vector x and returns the scalar function value, its gradient, and its Hessian,
% [functionValue,functionGradient,functionHessian] = f(x)
% the Matlab vectors x and lambda are the initial point and the initial Lagrange multiplier
% A and b are the matrix and vector of the linear constraint
% 1/2mu is the weight of the constraint penalisation in the augmented Lagrangian method
% maxNumSteps is maximum number of iterations
% display=true displays information about current iterates

function [x,lambda] = augmentedLagrangianMethod( f, x, lambda, A, b, mu, maxNumSteps, display, displayfun )
  % specify options of inner iterations
  if display
    displayval = 'iter';
  else
    displayval = 'off';
  end
  options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','MaxIter',30,'OptimalityTolerance',1e-5,'FunctionTolerance',1e-12,'Display',displayval);
  %options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxIter',1000,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-12,'Display',displayval);
  
  % optimization 
  iter = 0;
  constraintViolation = 1;
  tol = 1e-5;
  while ( iter <= maxNumSteps && constraintViolation > tol )
    
    % solve for x
    fun = @(x) augmentedLagrangian(f,A,b,mu,lambda,x);
    [x,~] = fminunc( fun, x, options );
    %x = newtonMethod( fun, x, 10, 1e-9, true, @(x) displayfun(x,lambda) );
    
    
    % update lambda and mu
    constraint = A*x-b;
    lambda = lambda - constraint/mu;
    mu = 0.99*mu;
    
    % output
    constraintViolation = norm(constraint,2);
    iter = iter+1;
    if display
      fprintf('finished outer iteration %d, constraint violation %e\n',iter,constraintViolation);
      displayfun(x,iter);
    end
  end
end

% augmented Lagragian: f(x)-lambda'(Ax-b)+1/(2*mu)||Ax-b||^2
function [e,de,d2e] = augmentedLagrangian( f, A, b, mu, lambda, x )
  constraint = A*x-b;
  [e,de,d2e] = f(x);
  e = e-lambda'*constraint+constraint'*constraint/(2*mu);
  de = de+A'*(-lambda+constraint/mu);
  d2e = d2e+1/mu*A'*A;
end
