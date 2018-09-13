function main(n,smooth,example,val) 
    %%%%%% INIZIALIZATION %%%%%%
    % Scelgo tau in funzione di tre parametri a0,a1,b1 e f deriva direttamente
    % dalle scelte (Per ottenere questa f e' necessario scegliere K=b1)
    % tau = @(x) min(a(i)*x+b(i))
    epsilon = .2;
    delta = 1e-5; % parameter that is used to approximate the norm i.e. we consider sqrt(sigma1^2 +sigma^2 + delta)...

    global gIter;
    gIter = 1;
    
    % Mesh initialization
    x = linspace(0,1,n);                                                    % Domain
    dx = x(2)-x(1);                                                         % Mesh size
    [X,Y] = meshgrid(x,x);
    Tri = delaunayTriangulation( X(:),Y(:) );                               % Triangularization
    nodes = Tri.Points;                                                     % Nodes of the mesh
    elements = Tri.ConnectivityList;                                        % Elements of the mesh
    bd_elements = Tri.freeBoundary;                                         % Constraint elements
    in_elements = setdiff([1:n^2]',bd_elements(:,1));                       % Interior elements

   % Utility matrices
    [M,S,H] = massStiffSquaresMatrix(nodes,elements);                       % Mass, Stiff and Squares matrices associated to the mesh
    S_in = S(in_elements,in_elements);                                      % Stiff matrix reduced to interior points
    M_in = M(in_elements,in_elements);                                      % Mass matrix reduced to interior points
    H_in = H(in_elements,in_elements);                                      % Squares matrix reduced to interior points
    nIn_el = numel(in_elements);                                            % Number of internal elements
    nBd_el = n^2 - nIn_el;                                                  % Number of boundary elements      
    
    % Matrices divergence constraint
    M1 = mixedMassStiffMatrix(nodes,elements,1);
    M2 = mixedMassStiffMatrix(nodes,elements,2);
    Aeq = [ M1' M2'];                                                       % linear divergence constraint for sigma
    
    % Constraint initialization
    [beq,diffMeasure,a,b] = constraint(n,x,example,dx,smooth,M);
    f_fun = @(x) f_function(x,a,b);
    
    dir = saveData(val, n,a,b,epsilon,delta,example,smooth);              % Routine to save the data

    img1 = figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(flip(diffMeasure)); colorbar; axis square;                      % Constraint visualization
    saveas(img1,[dir, '/diffMeasure.png'])
    close(img1);
    
    
    % Initial guess for phi and sigma 
    phi0 = .5*(b(1)+b(end))*ones(n^2,1); % rand(2*n^2,1) 
    sigma0 = zeros(2*n^2,1);
    
    
    %%%%%% OPTIMIZATION %%%%%%
    % simulated annealing iteration: decrease epsilon stepwise
    lambda = zeros(size(beq));
    while ( epsilon > dx )
        fprintf('Value of epsilon: %f\n',epsilon);
        % Intialize functions
        phiSolver = @(sigma,phi) PHIsolver(sigma,phi0,f_fun,b(1),M_in,S_in,H_in,delta,epsilon,in_elements,n);
        sigmaEnergy = @(sigma,phi) discreteEnergy( nodes, elements, phi,sigma, f_fun, delta );
        
        [sigma1,phi1,lambda] = minProcedure(phi0,sigma0,lambda,sigmaEnergy,phiSolver,Aeq,beq,.0001, 5, true,@(sigma,phi) visualize(sigma,phi,x,n,epsilon,dir,f_fun,b(1)));
        phi0 = phi1;
        sigma0 = sigma1;
        
        epsilon = epsilon / 2;                                              % Decrease epsilon

    end
    
    
    
    
end


%%%%%%  ROUTINE TO SAVE ALL DATA O  F THE SIMULATION %%%%%%
function dir = saveData(val, n,a,b,epsilon,delta,example,smooth)
    dir = ['Simulations/Results_',date,'_',num2str(val)];
    mkdir(dir);
    F = fopen([dir,'/data.txt'],'w');
    fprintf(F,['n = ',num2str(n) ,'\n']);
    fprintf(F,['a = ',num2str(a') ,'\n']);
    fprintf(F,['b = ',num2str(b') ,'\n']);
    fprintf(F,['epsilon = ',num2str(epsilon) ,'\n']);
    fprintf(F,['delta = ',num2str(delta) ,'\n']);
    fprintf(F,['example = ',num2str(example) ,'\n']);
    fprintf(F,['smooth = ',char(string(smooth))]); % ,' & radius dsmt = ', num2str(dsmt) ,'\n']);
    fclose(F);
end


%%%%%% VISUALIZATION FUNCTION %%%%%%
function visualize(sigma,phi,x,n,epsilon,dir,f,b0)
    global gIter;
    sigma = reshape(sigma(:),n,n,2);
    sigmaNorm = sqrt(sigma(:,:,1).^2+sigma(:,:,2).^2);
    phi = reshape(phi(:),n,n);
    [fphi,~] = f(phi(:));
    fphi= reshape(fphi,n,n);
    
    img2 = figure('units','normalized','outerposition',[0 0 1 1]);%,'pos',[50,80,1600,800]);
	quiver(x,x,sigma(:,:,1),sigma(:,:,2)); axis image;
    title(['\epsilon = ',num2str(epsilon) ]);
    saveas(img2,[dir, '/phi_sigma_at_iter', num2str(gIter),'.png']);
    close(img2); 

    img3 = figure('units','normalized','outerposition',[0 0 1 1]);%,'pos',[50,80,1600,800]);
    subplot(2,2,1); contourf(x,x,fphi); colorbar; axis image;     title(['\epsilon = ',num2str(epsilon) ]);    
    subplot(2,2,2); contourf(x,x,phi); colorbar; axis image;
    subplot(2,2,3); contourf(x,x,sigmaNorm); colorbar; axis image; 
    subplot(2,2,4); contourf(x,x,fphi.*sigmaNorm+(b0-phi).^2); colorbar; axis image;
    saveas(img3,[dir, '/fphi_phi_at_iter', num2str(gIter),'.png']);
    close(img3); 
    gIter = gIter+1;
end
