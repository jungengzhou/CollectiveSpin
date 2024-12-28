 

classdef DickeTools < handle  
    properties    
        N  
        Jx
        Jy
        Jz
        Jp
        Jn
        Ja
    end  
      

    methods  

        function obj = DickeTools(N)  
            obj.N = N;  
           [obj.Jx,obj.Jy,obj.Jz] = J_operator(N);
            obj.Jp = obj.Jx + 1i* obj.Jy;
            obj.Jn = obj.Jx  - 1i* obj.Jy;

            import containers.Map;  
            Ja = containers.Map;  
            Ja("x") = obj.Jx;
            Ja("y") = obj.Jy;
            Ja("z") =  obj.Jz;

            obj.Ja = Ja;
        end 

        function U = Rx(obj, phi)  
            U =  expm(-1i*phi*obj.Jx);
        end

        function U = Ry(obj, phi)  
            U =  expm(-1i*phi*obj.Jy);
        end

        function U = Rz(obj, phi)  
            U =  expm(-1i*phi*obj.Jz);
        end

        function U = Ra(obj, phi,direction)  
            if nargin == 1
                direction = "z";
            end
            U =  expm(-1i*phi*obj.Ja(direction));
        end

        function U = OAT(obj, chit,direction)  
            if nargin == 1
                direction = "z";
            end
            U =  expm(-1i*chit*obj.Ja(direction)^2);
        end

        function U = TNT(obj, chit,omega,direction)  
            if nargin == 1
                direction = "z";
            end
            U =  expm(-1i*chit* (obj.Ja(direction)^2+ omega*obj.Jx));
        end


        function state = SCS(obj,theta,phi)
           state = zeros(obj.N+1,1);
           state(1,1) = 1; 
           state = obj.Rz(phi) * obj.Ry(theta)*state;
        end

        function state = TFS(obj)
           state = zeros(obj.N+1,1);

           if  mod(obj.N, 2) == 0
               state(obj.N/2+1,1) = 1; 
           elseif mod(obj.N, 2) == 1
               ind= floor(obj.N/2 ) + 1;
                state(ind,1) = 0.5^0.5;
                state(ind+1,1) = 0.5^0.5;
           end
        end

        function state = CAT(obj,theta,phi)
           state = obj.SCS(theta,phi) + obj.SCS(pi-theta,phi);
           n = sum( abs(state).^2 );
           state = state/(n^0.5);
        end

        function Ob = Parityb(obj,state)
                        
            Ns = linspace(0,obj.N,obj.N+1)';
            Nbs = exp(1i*pi*Ns);
             
            O = abs(state).^2.*Nbs;
          
            Ob = real(sum(O));
        end

        % function prob = state_distribution(obj,state):
        % 


   end


end



function [Jx,Jy,Jz] = J_operator(N)

% 模a表示自旋上，模b表示自旋下
Dimentions = N+1;
adegab_vector = zeros(N,1);
bdegaa_vector = zeros(N,1);
adegaa_vector = zeros(N+1,1);
bdegab_vector = zeros(N+1,1);
%%
for Na=Dimentions-1:-1:1
    bdegaa_vector(N-Na+1) = sqrt(Na)*sqrt(N-Na+1);%Na从N到1
    adegab_vector(N-Na+1) = sqrt(Na)*sqrt(N-Na+1);%Na从N-1到0
end
%%
for ii=1:Dimentions
    adegaa_vector(ii) = N+1-ii;
    bdegab_vector(ii) = ii-1;
end
%%
bdegaa = diag(bdegaa_vector,-1);
adegab = diag(adegab_vector,+1);
adegaa = diag(adegaa_vector);
bdegab = diag(bdegab_vector);

%Jx,Jy,Jz的表示
Jx = (1/2)*(bdegaa + adegab);
Jy = (1i/2)*(bdegaa - adegab);
Jz = (1/2)*(adegaa - bdegab);
end



