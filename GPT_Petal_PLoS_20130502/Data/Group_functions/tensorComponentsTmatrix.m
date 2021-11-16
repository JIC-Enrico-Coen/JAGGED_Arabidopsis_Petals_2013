function [Glocal,J, Gmatrix] = tensorComponentsTmatrix( Gglobal )
%[Glocal,J] = tensorComponents( Gglobal )
%   Gglobal is a growth tensor represented as a 3-vector or a 6-vector.
%   Glocal is set to the three principal components of growth in descending
%   order and J to the corresponding orthonormal frame.  The columns of J
%   are the principal axes of the tensor.

    
    symmetrycount = 2;
    if length(Gglobal)==3
        Gmatrix = [ Gglobal(1) Gglobal(3)/symmetrycount; ...
                    Gglobal(3)/symmetrycount Gglobal(2) ];
    else
        Gmatrix = [ Gglobal(1) Gglobal(6)/symmetrycount Gglobal(5)/symmetrycount; ...
                    Gglobal(6)/symmetrycount Gglobal(2) Gglobal(4)/symmetrycount; ...
                    Gglobal(5)/symmetrycount Gglobal(4)/symmetrycount Gglobal(3) ];
    end

    [J,D] = eig( Gmatrix );
    Glocal = diag(D)';

    if length(Gglobal)==3
        Glocal = Glocal([2 1]);
        J = J(:,[2 1]);
    else
        Glocal = Glocal([3 2 1]);
        J = J(:,[3 2 1]);
    end
  % check = Gmatrix*J - J*diag(Glocal)
end
