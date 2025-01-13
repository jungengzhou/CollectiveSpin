function Q = HusimiQ_spin(state)
    % The HusimiQ distribution to visualize the states in collective spin 
    % By Jungeng Zhou, updated 2025-01-13.
    % https://github.com/jungengzhou/CollectiveSpin/tree/main
    
    shape = size(state);
    N = shape(1,1)-1;
    sys = CollectiveSpin(N);
    
    n = 100;     % grid points
    theta = linspace(0,pi,n);
    phi = linspace(0,2*pi,n);
    SCS = zeros(N+1,1,n,n);
    
    Q = zeros(n,n);
    wb=waitbar(0,'please wait');
    for i = 1:n
        for j = 1:n
            SCS(:,:,i,j) = sys.SCS(theta(i),phi(j));
            rhoc = SCS(:,:,i,j)*SCS(:,:,i,j)';
            Q(i,j) = trace(rhoc*state); % Husimi-Q function for density matrix
        end
        str=['calculation...',num2str(i/n*100),'%'];
        waitbar(i/n,wb,str)
    end
    delete(wb);
    Q = real(Q);
end

