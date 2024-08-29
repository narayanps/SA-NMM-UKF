function [y_new] = nmm_jr_ukf_C2C4(y, p)
% ------------------------------------------------------------------------------
% Author: Narayan P. Subramaniyam
% Affiliation: MET Faculty, Tampere University
% Email: [narayan.subramaniyam@tuni.fi]
%
% Description:
% This MATLAB code is developed as part of my research. Feel free to reuse 
% or modify this code, provided that you give proper attribution by citing 
% the associated paper. 
% 
%
%
% License:
% This code is licensed under a Creative Commons Attribution 4.0 International License.
% You are free to share and adapt the material for any purpose, even commercially,
% under the following terms:
% 1. You must give appropriate credit, provide a link to the license, and indicate 
%    if changes were made.
% 2. You must cite the original paper if you use this code in your work.
%
% For more details on the license, visit:
% https://creativecommons.org/licenses/by/4.0/
% ------------------------------------------------------------------------------

a=100;
A=6;
B=70;
b=50;
dt=0.001;
C(1) = 135;
C(2)=y(7);
C(3) = 33.75;
C(4)=y(8);

v_0=6;
v_P=v_0;
v_I=v_0;
G=0;
%pyramidal
y_new(1) = euler_f1(y(1), y(4), dt);
D = A*a*S(y(2)-y(3),v_P) ;
y_new(4) = euler_f2(D, a, y(4), y(1), dt );



%secondary pyramidal
y_new(2) = euler_f1(y(2), y(5), dt);
v_P_prime=v_0;
D = A*a*C(2)*S(C(1)*y(1),v_P_prime) + A*a*G*S(y(2)-y(3), v_P) + A*a*p;
y_new(5) = euler_f2(D, a, y(5), y(2), dt );

%inhibitory
y_new(3) = euler_f1(y(3), y(6), dt);
D = B*b*C(4)*S(C(3)*y(1),v_I) ;
y_new(6) = euler_f2(D, b, y(6), y(3), dt );
y_new(7) = C(2);
y_new(8) = C(4);
y_new=y_new';

function [z] = S(x,v)
        e0=2.5;
        r=0.56;
        z=2*e0/(1+exp(r*(v-x)));
end


function [y] = euler_f1(y,x,dt)
        y = y + dt*x;
end


function [y] = euler_f2(K,c,y,z,dt)
        y = y + dt*(K - 2*c*y - c^2*z);
end



function [z] = H(x,k)
        z=x/(x+k);
end

end

