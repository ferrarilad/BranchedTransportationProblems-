function [sigma,phi,lambda] = minProcedure(phi0,sigma0,lambda,energy,phi_solver,A,b,mu, maxNumSteps, display, displayfun )
    global gIter;
    % specify options of inner iterations
    if display
    displayval = 'iter';
    else
    displayval = 'off';
    end
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','MaxIter',20,'OptimalityTolerance',1e-5,'FunctionTolerance',1e-12,'Display',displayval);
    %options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxIter',100,'OptimalityTolerance',1e-5,'FunctionTolerance',1e-12,'Display',displayval);

    % optimization 
    iter = 0;
    constraintViolation = 1;
    tol = 1e-5;
    
    % initialize
    phi = phi0;
    while ( iter <= maxNumSteps && constraintViolation > tol )
        fprintf(['Global iteration number: ' , num2str(gIter)]);
        % solve for x
        f = @(sigma) energy(sigma,phi);
        fun = @(sigma) augmentedLagrangian(f,A,b,mu,lambda,sigma);
        
        [sigma,~] = fminunc( fun, sigma0, options );  
        phi = phi0;
        for j=1:5
            phi = phi_solver(sigma,phi);
        end
        % update lambda and mu
        constraint = A*sigma-b;
        lambda = lambda - constraint/mu;
        mu = 0.99*mu;

        % output
        constraintViolation = norm(constraint,2);
        iter = iter+1;
        if display
          fprintf('finished outer iteration %d, constraint violation %e\n',iter,constraintViolation);
          displayfun(sigma,phi);
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
