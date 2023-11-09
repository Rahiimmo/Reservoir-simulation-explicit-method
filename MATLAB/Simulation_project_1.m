clear
clc
close all
miu = 1        ;                                                           % viscosity              [cp]
c = 10 ^(-6)   ;                                                           % copressbility          [1/psi]
k = 100        ;                                                           % permeability           [md]
phi = 0.18     ;                                                           % porosity
P_0 = 0        ;                                                           % initial pressure         [psi]
L = 1          ;                                                           % total length             [m]
P_L = 100      ;                                                           % Left hand side pressure  [psi]
P_R = 0        ;                                                           % Right hand side pressure [psi]

% field units to SI units
miu = miu / 1000  ;                                                        % unit = pa.s
c = c * 14.7 / 101325   ;                                                  % unit = 1/pa
k = k * 9.869233 * 10^(-16) ;                                              % unit = m^2
P_L = P_L * 101325 / 14.7   ;                                              % unit = pa
alpha = (miu * c * phi) / k   ;                                            % unit = s/m^2

% assume:  Δx = dx   , Δt = dt
 delta_x = [ 0.01 , 0.1 , 0.2 ]  ;
 delta_t = [ 10^(-5) , 10^(-4) , 10^(-4) ]  ;
 T1 = []   ;
 T2 = []   ;
 T3 = []   ;
 T_list = { T1 , T2 , T3 }   ;
 X1 = []   ;
 X2 = []   ;
 X3 = []   ;
 X_list = { X1 , X2 , X3 }   ;
 S = [ -1 , -1 , -1 ]        ;
 
%% explicit method

P_1_ex = {}   ;
P_2_ex = {}   ;
P_3_ex = {}   ;
P_list_ex = { P_1_ex , P_2_ex , P_3_ex }   ;

for i = ( 1 : length(delta_x) )
   dx = delta_x(i)   ;
   dt = delta_t(i)   ;
   N = L / dx   ;
   T = T_list{i}   ;
   X = X_list{i}   ;
   P = P_list_ex{i}   ;
   s = 0.1 / dt   ;
   S(i) = s   ;
   beta = alpha / dt * dx^2   ;
   
   for j = ( 0 : s )
       tt = j * dt  ;
       T = [ T , tt ]   ;
       P{i+1} = []   ;
   end
   
   for j = (0:N-1)
       xx = (j + 0.5 ) * dx   ;
       X = [ X , xx ]   ;
       
       for r = ( 1 : s+1)
           P{r}(j+1) = 0   ;
       end
   end
   
   X = [ 0 , X , 1 ]   ;
   
   for j = ( 2 : s + 1 ) 
       
       P{j}(1) = P{j-1}(1) + 4 / ( 3 * beta ) * ( P{j-1}(2) - 3 * P{j-1}(1) + 2 * P_L )   ;
       
       for z = ( 2 : N - 1 )
           P{j}(z) = P{j-1}(z) + 1 / beta * ( P{j-1}(z+1) - 2 * P{j-1}(z) + P{j-1}(z-1) )  ;
       end
       
       P{j}(N) = P{j-1}(N) + 4 / ( 3 * beta ) * ( 2 * P_R - 3 * P{j-1}(N) + P{j-1}(N-1) )   ;
          
   end
   
   for z = ( 1 : s + 1 ) 
       P{z} = [ P_L , P{z} , P_R ]   ;
   end
   
   T_list{i} = T   ;
   X_list{i} = X   ;
   P_list_ex{i} = P  ;
    
end

for i = ( 1 : length(delta_x) )
    for j = ( 0 : log10( S(i) - 1 ) ) 
        for q = ( 1 : 9 )
            z = q * 10^j  ;
            figure(i)
            title(['explicit - Δx = ',num2str( 100 * delta_x(i) ),'cm'])
            hold on
            plot ( X_list{i} , P_list_ex{i}{z} )  ;
            ylabel('Pressure(pa)')  ;
            xlabel('x(m)')  ;
        end
    end
end

%% Implicit method

P_1_im = {}   ;
P_2_im = {}   ;
P_3_im = {}   ;
P_list_im = { P_1_im , P_2_im , P_3_im }   ;

for i = ( 1 : length(delta_x) )
   dx = delta_x(i)   ;
   dt = delta_t(i)   ;
   N = L / dx   ;
   P = P_list_im{i}   ;
   s = S(i)   ;
   beta = alpha / dt * dx^2   ;
   
   for j = ( 1 : s + 1 )
       P{j} = []   ;
       for j1 = ( 1 : N )
           P{j}(j1) = 0  ;
       end
   end
   
  A = zeros (N)   ;
  B = zeros ([ N , 1 ])   ;
   
   for j = ( 2 : s + 1 ) 
       A(1,1) = - 3 / 4 * alpha - 3   ;
       A(1,2) = 1   ;
       B(1) = - 3 / 4 * alpha * P{j-1}(1) - 2 * P_L   ;
       
       for z = ( 2 : N - 1 )
           A(z,z-1) = 1   ;
           A(z,z+1) = 1   ;
           A(z,z) = - alpha - 2   ;
           B(z) = - alpha * P{j-1}(z)   ;
       end
       
       A(N,N-1) = 1   ;
       A(N,N) = - 3 / 4 * alpha - 3   ;
       B(N) = - 3 / 4 * alpha * P{j-1}(N) - 2 * P_R   ;
       
       P{j} = (linsolve(A,B))'   ;  
   end
   
   for z = ( 1 : s + 1 ) 
       P{z} = [ P_L , P{z} , P_R ]   ;
   end
   
   P_list_im{i} = P  ;
  
end

for i = ( 1 : length(delta_x) )
    for j = ( 0 : log10( S(i) - 1 ) ) 
        for q = ( 1 : 9 )
            z = q * 10^j  ;
            figure(i+ length(delta_x))
            title(['implicit - Δx = ',num2str( 100 * delta_x(i) ),'cm'])
            hold on
            plot ( X_list{i} , P_list_im{i}{z} )  ;
            ylabel('Pressure(pa)')  ;
            xlabel('x(m)')  ;
        end
    end
end

%% Crank-Nicholson method

P_1_CN = {}   ;
P_2_CN = {}   ;
P_3_CN = {}   ;
P_list_CN = { P_1_CN , P_2_CN , P_3_CN }   ;

for i = ( 1 : length(delta_x) )
   dx = delta_x(i)   ;
   dt = delta_t(i)   ;
   N = L / dx   ;
   P = P_list_CN{i}   ;
   s = S(i)   ;
   beta = alpha / dt * dx^2   ;
   
   for j = ( 1 : s + 1 )
       P{j} = []   ;
       for j1 = ( 1 : N )
           P{j}(j1) = 0  ;
       end
   end
   
  A = zeros (N)   ;
  B = zeros ([ N , 1 ])   ;
   
   for j = ( 2 : s + 1 ) 
       A(1,1) = - 3 / 2 * beta - 3   ;
       A(1,2) = 1   ;
       B(1) = ( - 3 / 2 * beta + 3 ) * P{j-1}(1) - P{j-1}(2) - 4 * P_L   ;
       
       for z = ( 2 : N - 1 )
           A(z,z-1) = 1   ;
           A(z,z+1) = 1   ;
           A(z,z) = - 2 * beta - 2   ;
           B(z) = - P{j-1}(z-1) + ( - 2 * beta + 2 ) * P{j-1}(z) - P{j-1}(z+1)   ;
       end
       
       A(N,N-1) = 1   ;
       A(N,N) = - 3 / 2 * beta - 3   ;
       B(N) = - P{j-1}(N-1) + ( - 3 / 2 * beta + 3 ) * P{j-1}(N) - 4 * P_R   ;
       
       P{j} = (linsolve(A,B))'   ;  
   end
   
   for z = ( 1 : s + 1 ) 
       P{z} = [ P_L , P{z} , P_R ]   ;
   end
   
   P_list_CN{i} = P  ;
  
end

for i = ( 1 : length(delta_x) )
    for j = ( 0 : log10( S(i) - 1 ) ) 
        for q = ( 1 : 9 )
            z = q * 10^j  ;
            figure(i+ 2 * length(delta_x))
            title(['Crank-Nicholson - Δx = ',num2str( 100 * delta_x(i) ),'cm'])
            hold on
            plot ( X_list{i} , P_list_CN{i}{z} )  ;
            ylabel('Pressure(pa)')  ;
            xlabel('x(m)')  ;
        end
    end
end
 
%% Analytical method

P_1_A = {}   ;
P_2_A = {}   ;
P_3_A = {}   ;
P_list_A = { P_1_A , P_2_A , P_3_A }   ;

for i = ( 1 : length(delta_x) )
   dx = delta_x(i)   ;
   dt = delta_t(i)   ;
   N = L / dx   ;
   T = T_list{i}   ;
   X = X_list{i}   ;
   P = P_list_A{i}   ;
   s = S(i)  ;
   
   for j = ( 1 : s + 1 )
       P{j} = []   ;
       for j1 = ( 1 : N )
           P{j}(j1) = 0  ;
       end
   end
   
   for j = ( 2 : s + 1 ) 
       for z = ( 2 : N + 1 )
           x = X(z)   ;
           t = T(j)   ;
           n = 1   ;
           sum = 0   ;
           while ( ( 1 / n ) * exp( - t * n^2 * pi^2 / ( L^2 * alpha ) ) * sin( n * pi * x / L ) ) > 0.000001
               sum = sum + ( 1 / n ) * exp( - t * n^2 * pi^2 / ( L^2 * alpha ) ) * sin( n * pi * x / L )   ;
               n = n + 1    ; 
           end
           
           P{j}(z-1) = P_L + ( P_R - P_L ) * ( x / L + 2 / pi * sum )   ;
           if P{j}(z-1)<0
               P{j}(z-1) = 0   ;
           end
      
       end
   end
   
   for z = ( 1 : s + 1 ) 
       P{z} = [ P_L , P{z} , P_R ]   ;
   end
   
   P_list_A{i} = P   ;
   
end
 
for i = ( 1 : length(delta_x) )
    for j = ( 0 : log10( S(i) - 1 ) ) 
        for q = ( 1 : 9 )
            z = q * 10^j  ;
            figure(i+ 3 * length(delta_x))
            title(['Analytical - Δx = ',num2str( 100 * delta_x(i) ),'cm'])
            hold on
            plot ( X_list{i} , P_list_A{i}{z} )  ;
            ylabel('Pressure(pa)')  ;
            xlabel('x(m)')  ;
        end
    end
end
 
 
 
 
 

