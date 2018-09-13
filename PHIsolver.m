% Function to solve the PDE in weak form associated to the optimization problem
% int f(phi)|sigma| + 1/2[epsilon |grad phi|^2+(phi-b(1))^2/epsilon] dx
% with Dirichelet boundary data equals to b(1).  
function phi = PHIsolver(sigma,phi0,f_fun,b,M_in,S_in,H_in,delta,epsilon,in_elements,n)
    % f is piecewise affine so f_der is piecewise affine. Divide on each region and then the problem is solvable 
    [~,f_der] = f_fun(phi0(in_elements));                                  % Evaluation of the gradient on the previous optimal region
    sigma1 = sigma(1:n^2);
    sigma2 = sigma(n^2+1:end);
    regSigmaNrm = sqrt(sigma1(in_elements).^2 + sigma2(in_elements).^2+delta);
    phi = zeros(n^2,1);
    A = epsilon*S_in+1/epsilon*M_in;
    B = -H_in*f_der.*regSigmaNrm;
    
    % Optimize --> int f(phi)|sigma| + 1/2[epsilon*|graphi|^2+(phi-b0)^2/epsilon] dx
    phi_in = A\B;                                                           % by exactly solving the associated PDE in weak form.
    
    % Assemble solution w/ boudary constraint
    phi(in_elements) = phi_in; %max(min(phi_in+1,1),0);
    phi = phi+b; %(bd_elements(:,1)) = ones(nBd_el,1); 

end
