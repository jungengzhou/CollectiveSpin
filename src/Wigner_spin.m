function W = Wigner_spin(Op)
    % The Wigner distribution of collective spin  (j=N/2) operator -- Op
    % To obtain the Wigner distribution, the operator is expanded in terms of spherical tensors
    % T_{k,q} = \sum_{m,m'=-j}^{j}{ (-1)^(j-m)*(2*k+1)^(1/2)*Wigner3j([j,k,j],[-m,q,m'] ) |m><m'|   }
    % c_{k,q} = Trace(Op*T_{k,q})
    % W(theta,phi) = \sum_{k=0}^{N}{  \sum_{q=-k}^{k} { c_{k,q}*harmonicY_(k,q,theta,phi) } }

    % Wigner3j is Wigner 3j symbol, which is referenced from Kobi. 
    % Wigner3j symbol (https://www.mathworks.com/matlabcentral/fileexchange/20619-wigner3j-symbol), MATLAB Central File Exchange
    % harmonicY_(k,q,theta,phi) is the Spherical harmonic function, which is referenced from Javier Montalt Tordera
    % Spherical Harmonics (https://github.com/jmontalt/harmonicY/releases/tag/v2.0.1), GitHub

    % By Jungeng Zhou, 2025-01-13.
    % https://github.com/jungengzhou/CollectiveSpin/tree/main

    % Get the shape of the input Op or operator (matrix)
    shape = size(Op);
    % Determine the article number based on the matrix size
    N = shape(1,1)-1;

    % Define the number of grid points for theta and phi
    n = 100;     
    thetas = linspace(0,pi,n);
    phis = linspace(0,2*pi,n);

    % Initialize the Wigner function matrix
    W = zeros(n,n);

    % Calculate all ckq coefficients for the given Op
    ckqs = ckq_all(Op);

    % Create a waitbar to show progress
    wb=waitbar(0,'please wait');
    for i = 1:n
        for j = 1:n
            % Retrieve theta and phi values for the current grid point
            theta = thetas(i);
            phi = phis(j);
            % Calculate the spherical harmonics for all m values
            M = ckqs.*harmonicY_all(N,theta,phi);
            % Sum the product of ckqs and spherical harmonics for the current grid point
            W(i,j) = W(i,j) + sum(M(:));
        end
        % Update the waitbar with the current progress
        str=['calculation...',num2str(i/n*100),'%'];
        waitbar(i/n,wb,str)
    end
    % Delete the waitbar after the loop is complete
    delete(wb);

    % Convert the Wigner function matrix to real values
    W = real(W);
end

function harmonicYs = harmonicY_all(N,theta,phi)
    % Initialize the matrix to store spherical harmonics
    harmonicYs = zeros(N+1,2*N+1);
    
    % Loop through all possible k and q values
    for k = 0:N
        for q = -k:k
            % Calculate the spherical harmonic Y(k,q,theta,phi)
            harmonicYs(k+1,k+q+1) = harmonicY(k,q, theta, phi);
        end
    end

end

function ckqs = ckq_all(Op)
    % Get the shape of the input Op matrix
    shape = size(Op);
    % Determine the maximum quantum number based on the matrix size
    N = shape(1,1)-1;
    % Initialize the matrix to store ckq coefficients
    ckqs = zeros(N+1,2*N+1);

    % Loop through all possible k and q values
    for k = 0:N
        for q = -k:k
            % Calculate the ckq coefficient for the given k and q
            ckqs(k+1,k+q+1) = ckq_matrix(Op,k,q);
        end
    end

end

function ckq = ckq_matrix(Op,k,q)
    % Get the shape of the input Op matrix
    shape = size(Op);
    % Determine the maximum quantum number based on the matrix size
    N = shape(1,1)-1;
    % Generate a linspace for quantum numbers m
    Ns = linspace(N/2,-N/2,N+1);

    % Define quantum numbers for the Wigner 3j symbols
    J1 = N/2; 
    J2 = k; M2 = q;
    J3 = N/2; 

    % Initialize the Tkq matrix
    Tkq = zeros(N+1,N+1);

    % Loop through all possible m and m_star values
    for i = 1:N+1
        for j = 1:N+1
            m = Ns(i);
            m_star = Ns(j);
            % Combine quantum numbers for the Wigner 3j symbols
            J123 = [J1, J2, J3];
            M123 = [-m, M2, m_star];       

            % Calculate the Wigner 3j symbol
            Tkq(i,j) = (-1)^(N/2-m)*(2*k+1)^0.5*Wigner3j( J123, M123 );
        end
    end

    % Calculate the ckq coefficient by taking the trace of Tkq*Op
    ckq = trace(Tkq*Op);

end