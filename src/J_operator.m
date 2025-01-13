
function [Jx,Jy,Jz] = J_operator(N)
% Three basic Collective spin operators Jx,Jy,Jz
% Inspired by Hongtao Huo
% By Jungeng Zhou, updated 2025-01-13.

Dimentions = N+1;

adaga_vector = linspace(N,0,N+1);
bdagb_vector = linspace(0,N,N+1);
adaga = diag(adaga_vector);
bdagb = diag(bdagb_vector);

adagb_vector = zeros(N,1);
for Na=Dimentions-1:-1:1
    adagb_vector(N-Na+1) = sqrt(Na)*sqrt(N-Na+1);
end

adagb = diag(adagb_vector,+1);
bdaga = adagb';
% bdaga = diag(adagb_vector,-1);

Jx = (1/2)*(bdaga + adagb);
Jy = (1i/2)*(bdaga - adagb);
Jz = (1/2)*(adaga - bdagb);

end

