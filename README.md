# Dry-coil-heat-exchanger-performance-evaluation-MATLAB-code-to-determine-HTC-and-fin-efficiency-
MATLAB code developed for determination of HTC and fin efficiency for a dry cross flow plain finned heat exchanger can be accomplished by using this code. Equivalent circular area of the plain fin, geometry inputs, Inlet &amp; outlet temperatures, heat transfer rate needs to be entered.  Solver technique: Succcessive Successive Substitution Method
%version: 1 dry coil d7_1P under moist air
%Author: MIDHUN V C (NYCU)
%Date: 24 July 2026
%Purpose: data reduction of dry coil
%Libraries used: 
%Reference Pirompugd et al 2005
function dry_coil ()
r1= 0.003795;   % metres (Tube outer radius)
r2= 0.012962925;   % metres (Equivalent fin radius / Schmidt r_eq)
ri= 0.0037;% pipe inner radius
k_f=236;% conduc of fin/collar
k_p=236;% 
t_fin=0.000095;%fin thickness
Ao=0.000876898; % total area in m2
A_pi=3.02221E-05;
A_po=2.87328E-05;% cal from outer radius and fin pitch
Lp=0.0013;% fin pitch or tube length
%A_pm=(A_po-A_pi)./log(A_po./A_pi);% mean area of the pipe for resist. term
A_f=0.000848166;
Xp=9.5E-05;%collar/pipe thickness ; used because of the mean area of the pipe A_pm%




T_a_in=20+273.15;
T_w=42.33+273.15;
Q=0.6586959;
T_a_out=33.048+273.15;
dT_in=abs(T_a_in-T_w);
dT_out=abs(T_a_out-T_w);
LMTD=(dT_in-dT_out)/(log(dT_in/dT_out));
R_total=LMTD/Q;
R_pipe=log(r1/ri)/(2*pi*k_p*Lp);
R_air=R_total-R_pipe;



%initialisation for the "for loop"


h_co=50;
del_hco= 1e-6;
iter_max=100;

% for Loop  %
for iter=1:iter_max


% expression for dry fin efficiency%
    
m=sqrt((2*h_co)/(k_f*t_fin));
mr1 = r1 * m; 
mr2 = r2 * m; 


    I0_r1 = besseli(0, mr1);
    I1_r1 = besseli(1, mr1);
    I1_r2 = besseli(1, mr2);
    K0_r1 = besselk(0, mr1);
    K1_r1 = besselk(1, mr1);
    K1_r2 = besselk(1, mr2);


    Nr=(K1_r1*I1_r2)-(K1_r2*I1_r1);
    Dr=(K1_r2*I0_r1)+(K0_r1*I1_r2);


    eta=((2*r1)/(m*(r2^2-r1^2)))*(Nr/Dr);
    
    
    h_co_new=1/(R_air*(A_po+(eta*A_f)));
    if abs(h_co-h_co_new)<del_hco
        fprintf('Converged at iteration=%d\n', iter);
        fprintf('Current h_co=%f\n',h_co);
        break
        
    end
    if iter ==100
    fprintf('Maximum iterations reached');
    end
h_co=h_co_new;% updating the solution
end

fprintf('Converged h_co=%.4f, eta = %.4f\n',h_co,eta);

end


