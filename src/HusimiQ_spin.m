function Q = HusimiQ_spin(state)
    % The HusimiQ distribution to visualize the states in collective spin 
    % By Jungeng Zhou, updated 2025-01-13.
    % https://github.com/jungengzhou/CollectiveSpin/tree/main
    
    shape = size(state);
    N = shape(1,1)-1;
    sys = CollectiveSpin(N);
    
    n = 100;     % grid points
    thetas = linspace(0,pi,n);
    phis = linspace(0,2*pi,n);
    
    Q = zeros(n,n);
    wb=waitbar(0,'please wait');
    for i = 1:n
        theta = thetas(i);
        parfor j = 1:n
            phi = phis(j);
            SCS = sys.SCS(theta,phi);
            rhoc = SCS*SCS';
            Q(i,j) = trace(rhoc*state); % Husimi-Q function for density matrix
        end
        str=['calculation...',num2str(i/n*100),'%'];
        waitbar(i/n,wb,str)
    end
    delete(wb);
    Q = real(Q);
end

