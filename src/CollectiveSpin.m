

classdef CollectiveSpin < handle  
    % Define the class of collective spin operators and states
    % By Jungeng Zhou, updated 2025-01-13. 
    % https://github.com/jungengzhou/CollectiveSpin/tree/main

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

        function obj = CollectiveSpin(N)  
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

        function phi=phi_OAT_opt(obj,chit)
        Bj = 4*(cos(chit))^(obj.N-2)*sin(chit);
        Aj = 1 - (cos(2*chit))^(obj.N-2);
        phi =  -0.5 * atan( Bj / Aj );
        end

         function U = TACT(obj, chit)  
            U =  expm(-1i*chit*(obj.Jy*obj.Jz+obj.Jz*obj.Jy) );
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

        function state = SSS_opt(obj,chit)
           state = obj.Rx(pi/2+obj.phi_OAT_opt(chit))*obj.OAT(chit,"z")*obj.SCS(pi/2,0); 
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




