function [dis, h] = bloch(state, distribution, views)
    % Visualization the state or operator with Husimi-Q or Winger distribution 
    % distribution: "Husimi"(default)  or "Wigner"
    % views: "X", "Y", "Z",  plot the 3D views when 
    % Thanks for the colormap provided by Yi Shen
    % By Jungeng Zhou, updated 2025-01-13.
    % https://github.com/jungengzhou/CollectiveSpin/tree/main

    % check the input type of state
    shape = size(state);
    if shape(2) == 1
        state = state*state';
    elseif shape(2) == shape(1)
    else
        error("Wrong input, the shape of state should be a vector or a square matrix");
    end

    % Set the default value of distribution to "Husimi"
    if nargin < 2 || isempty(distribution)
        distribution = "Husimi";
    end
    
    is2DView = false;

    % Determines if the views argument exists and is a character, and checks if it is one of 'X', 'Y', 'Z'
    if nargin > 2 && ischar(views)
        is2DView = any(strcmpi(views, {'X', 'Y', 'Z'}));
        if ~is2DView
            warning("Wrong views handling in code logic. Only X/Y/Z handling shown in example.");
        end
    end

    % Calculates the Husimi or Wigner distribution based on the distribution parameter
    if strcmpi(distribution, "Husimi")
        dis = HusimiQ_spin(state);
    elseif strcmpi(distribution, "Wigner")
        dis = Wigner_spin(state);
    else
        error("Wrong distribution in code, only Husimi/Wigner is supported.");
    end

    % Plot the Bloch Sphere
    h = myplot(dis, is2DView);
end


function h = myplot(Q, is2DView)
    n = 100;     % Husimi grid points
    theta = linspace(0, pi, n);
    phi = linspace(0, 2*pi, n);
    r = 1;  % Default radius for 3D plot
    
    figure;
    [phi2,theta2] = meshgrid(phi,-(theta-pi/2));
    [x,y,z] = sph2cart(phi2,theta2,r);
    
    % Plot the surface
    h = surf(x,y,z,Q, 'FaceAlpha', 0.6);
    shading interp;
    axis equal;
    grid on;
    axis off;
    hold on;
    
    gridcolor = '#8B8989';
    
    % Plot grid lines
    phiz = linspace(0,2*pi,n);
    for phix = linspace(0,2*pi,21)
        rz = r*cos(phiz);
        rx = r*sin(phiz)*cos(phix);
        ry = r*sin(phiz)*sin(phix);
        plot3(rx,ry,rz,Color=gridcolor,LineWidth=0.8,LineStyle="--");
    end
    
    phix = linspace(0,2*pi,n);
    for phiz = linspace(0,pi,21)
        rz = r*cos(phiz)*ones(1,n);
        rx = r*sin(phiz)*cos(phix);
        ry = r*sin(phiz)*sin(phix);
        plot3(rx,ry,rz,Color=gridcolor,LineWidth=0.8,LineStyle="--");
    end
    
    % Plot axes arrows
    L = r*1.5;
    if ~is2DView
        quiver3(0,0,0,1.7,0,0,'Color',"k",'LineWidth',2);
        quiver3(0,0,0,0,1.7,0,'Color',"k",'LineWidth',2);
        quiver3(0,0,0,0,0,1.6,'Color',"k",'LineWidth',2);
        view(135, 21);
        text(1.7, -0.2, 0.18, "Jx", 'FontSize', 20);
        text(-0.1, 1.4, 0.1, "Jy", 'FontSize', 20);
        text(0, 0.02, 1.58, "Jz", 'FontSize', 20);
    else
        quiver3(0,0,0,L,0,0,'Color',"k",'LineWidth',2);
        quiver3(0,0,0,0,L,0,'Color',"k",'LineWidth',2);
        quiver3(0,0,0,0,0,L,'Color',"k",'LineWidth',2);
        index = bin2dec(num2str(is2DView) );
        s = {'Z'; 'Y'; '' ;'X'};
        ViewStr = s(index);
        if strcmpi(ViewStr, 'X')
            view(90, 0);
            text(0, 1.25, 0.18, "Jy", 'FontSize', 25);
            text(0, 0.05, 1.4, "Jz", 'FontSize', 25);
        elseif strcmpi(ViewStr, 'Y') % Similarly, this block would be adjusted in a real scenario
            view(180, 0);
            text(1.4, 0, 0.14, "Jx", 'FontSize', 25);
            text(-0.03,0,  1.4, "Jz", 'FontSize', 25);   
        elseif strcmpi(ViewStr, 'Z')
            view(0, 90);
            text(1.2, 0.15, 0, "Jx", 'FontSize', 25);
            text(0.04, 1.4, 0, "Jy", 'FontSize', 25);
        else
            warning("Wrong views handling in code logic. Only X/Y/Z handling shown in example.");
        end
    end
    
    % Set colormap
    colormap(mymap());
end


function cmap = mymap()
    aaa = 42; % Number of colors in the colormap
    cmap = zeros(6*aaa, 3); % Create colormap matrix
    cmap(:, 3) = [linspace(1, 1, aaa) linspace(1, 1, aaa) linspace(1, 0, aaa) linspace(0, 0, aaa) linspace(0, 0, aaa) linspace(0, 0, aaa)]; % Blue
    cmap(:, 2) = [linspace(1, 0.75, aaa) linspace(0.75, 1, aaa) linspace(1, 1, aaa) linspace(1, 0.5, aaa) linspace(0.5, 0, aaa) linspace(0, 0, aaa)]; % Green
    cmap(:, 1) = [linspace(1, 0.83, aaa) linspace(0.83, 0.5, aaa) linspace(0.5, 1, aaa) linspace(1, 1, aaa) linspace(1, 1, aaa) linspace(1, 0, aaa)]; % Red
    colormap(cmap);
end