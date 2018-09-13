% nodes = n by 2 array whose rows are the mesh nodes
% elements = m by 3 array whose rows are the indices of the nodes of a triangle
% k is coordinate direction for derivative

function S = mixedMassStiffMatrix(nodes,elements,l)
% computes mixed mass and stiffness matrix ( \int_\Omega \phi_i \partial_l \phi_j dx )_{ij}
  n = length(nodes);
  rows = zeros(4*length(elements),1);
  cols = rows;
  vals = rows;

  I = eye(2);
  D = [-1 -1;I]';

  counter = 1;
  for k = 1:size(elements,1)          % loop through all elements {T}
      nodeInds = elements(k, 1:3);    % get global indices
      x = nodes(nodeInds, 1:2);       % local triangle vertices (as matrix rows)
      A = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
      AinvD = (D'/A)';                % gradient operator
      vol = abs(det(A)/2);            % triangle-volume

      % assemble matrices
      [c,r] = meshgrid(nodeInds,nodeInds);
      rows(counter:counter+8) = r(:);
      cols(counter:counter+8) = c(:);
      vals(counter:counter+8) = vol/3*[1 1 1]'*AinvD(l,:);
      counter = counter + 9;
  end

  S = sparse(rows,cols,vals,n,n);
end
