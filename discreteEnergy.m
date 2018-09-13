% computes the energy \int_\Omega f(\phi)|\sigma| dx and its first and second derivatives
% with respect to sigma. Uses linear finite elements with edge midpoint quadrature
function [e,de,d2e] = discreteEnergy( nodes, elements, phi,sigma, f_fun, delta )
    n = length(nodes);
    sigma = reshape(sigma(:),n,2);

    % set quadrature points on unit reference domain as well as Lagrange polynomial values
    % Lagrange polynomials on unit reference domain: l1=1-x-y, l2=x, l3=y
    quadPoints = [0 .5; .5 0; .5 .5]';
    quadWeights = [1 1 1]/3;
    valCoeffs = [1-quadPoints(1,:)-quadPoints(2,:);quadPoints(1,:);quadPoints(2,:)]';

    I = eye(2);
    D = [-1 -1;I]';

    % initialize vectors of row, column and value indices for second derivative matrix as well as e and de components
    rows = zeros(9*length(elements),1);
    cols = rows;
    vals1 = rows;
    vals2 = rows;
    vals3 = rows;

    e = 0;
    de = zeros(size(sigma));

    % precompute auxiliary variables
    auxMat1 = [];
    for j = 1:length(quadWeights)
        auxMat1 = [auxMat1, quadWeights(j) * kron(valCoeffs(j,:)',valCoeffs(j,:)')];
    end
    
    auxMat3 = bsxfun( @times, quadWeights', valCoeffs )';
    % assemble energy and derivatives
    counter = 1;
    for k = 1:size(elements,1)                                              % loop through all elements {T}
        nodeInds = elements(k, 1:3);                                        % get global node indices
        x = nodes(nodeInds, 1:2);                                           % triangle vertex positions (as matrix rows)
        A = (D * x)';                                                       % Jacobian of coordinate transform from T_ref to T.
        vol = abs(det(A)/2);                                                % triangle-volume
        sigmaVals = valCoeffs*sigma(nodeInds,:);                            % values of sigma at quadrature points
        phiVals = valCoeffs*phi(nodeInds);                                  % values of phi at quadrature points
        [fphi, ~] = f_fun(phiVals);                                         % values of f(phi) evaluated at the quadrature 
        regSigmaNrm = sqrt(sigmaVals(:,1).^2 + sigmaVals(:,2).^2+delta);    % value of the regularized norm of the vector measure sigma at quadrature points

        %%%%% provare %%%%%%%%%%%%%%%
        %[fphi1,~]=f_fun(phi(nodeInds));
        %fphi = valCoeffs*fphi1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % compute energy contribution on element
        e = e + vol * quadWeights * ( fphi .* regSigmaNrm);

        % compute derivative contribution on element
        de(nodeInds,1) = de(nodeInds,1) + vol * ( auxMat3 * ( fphi .* sigmaVals(:,1) ) ./ regSigmaNrm );
        de(nodeInds,2) = de(nodeInds,2) + vol * ( auxMat3 * ( fphi .* sigmaVals(:,2) ) ./ regSigmaNrm );

        % compute second derivative contribution on element
        [c,r] = meshgrid(nodeInds,nodeInds);
        rows(counter:counter+8) = r;
        cols(counter:counter+8) = c;

        vals1(counter:counter+8) = auxMat1 * ( vol * fphi .* ( sigmaVals(:,2).^2 + delta ) ./ regSigmaNrm.^3 );
        vals2(counter:counter+8) = auxMat1 * ( vol * fphi .* (-sigmaVals(:,1) .* sigmaVals(:,2) ) ./ regSigmaNrm.^3 );
        vals3(counter:counter+8) = auxMat1 * ( vol * fphi .* ( sigmaVals(:,1).^2 + delta ) ./ regSigmaNrm.^3 );

        counter = counter + 9;
    end

    % assemble second derivative matrix
    d2edsigma1dsigma1 = sparse(rows,cols,vals1,n,n);
    d2edsigma1dsigma2 = sparse(rows,cols,vals2,n,n);
    d2edsigma2dsigma2 = sparse(rows,cols,vals3,n,n);
    
    de = de(:);
    d2e = [  d2edsigma1dsigma1   d2edsigma1dsigma2  ;
             d2edsigma1dsigma2'  d2edsigma2dsigma2  ];
end


