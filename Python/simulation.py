import matplotlib.pyplot as plt
import scipy.linalg as la
import numpy as np
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']

miu = 1                                                                         # Viscosity                    [cp]
c = 10 ** (-6)                                                                  # Compressibility              [1 / psi]
k = 100                                                                         # Permeability                 [md]
phi = 0.18                                                                      # Porosity
P_0 = 0                                                                         # Initial pressure             [psi]
P_L = 100                                                                       # Left hand side pressure      [psi]
P_R = 0                                                                         # Right hand side pressure     [psi]
L = 1                                                                           # Total length                 [m]

# field units to SI units
miu = miu / 1000                                                                # unit = pa.s
c = c * 14.7 / 101325                                                           # unit = 1 / pa
k = k * 9.869233 * 10 ** (-16)                                                  # unit = m ** 2
alpha = (miu * c * phi) / k                                                     # unit = s / m ** 2
P_L = P_L * 101325 / 14.7                                                       # unit = pa

# assume: ∆x = dx  , ∆t = dt
delta_x = [ 0.01 , 0.1, 0.2 ]
delta_t = [ 10**(-5) , 10**(-4) , 10**(-4) ]
T1 = []
T2 = []
T3 = []
T_list = [ T1 , T2 , T3 ]
X1 = []
X2 = []
X3 = []
X_list = [ X1 , X2 , X3 ]
S = [ -1 , -1 , -1 ]

#explicit method

P_1_ex = []
P_2_ex = []
P_3_ex = []
P_list_ex = [ P_1_ex , P_2_ex , P_3_ex ]

for i  in range ( len ( delta_x ) ) :
    dx = delta_x[i]
    dt = delta_t[i]
    N = L / dx
    T = T_list[i]
    X = X_list[i]
    P = P_list_ex[i]
    s = 0.1 / dt
    S[i] = s

    for j in range ( int(s) + 1 ) :
        tt = j * dt
        T.append(tt)
        a = []
        P.append(a)

    for j in range ( int(N) ) :
        xx =  (j + 0.5) * dx
        X.append(xx)
        for r in range ( int(s) + 1 ) :
            P[r].append(0)

    X = [0] + X + [1]
    for j in range ( 1 , int(s) + 1 ) :
        P[j][0] = P[j-1][0] + ( 4 * dt ) / ( 3 * alpha * dx**2 ) * ( P[j-1][1] - 3 * P[j-1][0] + 2 * P_L )
        for z in range ( 1 , int(N-1) ) :
            P[j][z] = P[j-1][z] + dt / ( alpha * dx**2 ) * ( P[j-1][z+1] - 2 * P[j-1][z] + P[j-1][z-1] )
        P[j][int(N)-1] = P[j-1][int(N)-1] + ( 4 * dt ) / ( 3 * alpha * dx**2 ) * ( 2 * P_R - 3 * P[j-1][int(N)-1] + P[j-1][int(N)-2] )

    for z in range ( int(s) + 1 ) :
        P[z] = [P_L] + P[z] + [P_R]

    T_list[i] = T
    X_list[i] = X
    P_list_ex[i] = P

for i  in range ( len ( delta_x ) ) :
    for j in range ( int( np.log10(S[i]) ) ) :
        for q in [ 1 , 3 , 5 , 8 ] :
            z = q * 10**j
            plt.figure( 'explicit - ∆x = {}cm'.format( int( 100 * delta_x[i] ) ) )
            plt.title('explicit - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.xlabel('x(m)')
            plt.ylabel('Pressure(pa)')
            plt.plot( X_list[i] , P_list_ex[i][z] , label = '∆t = {}'.format( T_list[i][z] ) )
            plt.legend()


# implicit method
P_1_im = []
P_2_im = []
P_3_im = []
P_list_im = [ P_1_im , P_2_im , P_3_im ]

for i  in range ( len ( delta_x ) ) :
    dx = delta_x[i]
    dt = delta_t[i]
    N = L / dx
    P = P_list_im[i]
    s = S[i]

    for j in range ( int(s) + 1 ) :
        a = []
        P.append(a)
        for j1 in range ( int(N) ) :
            P[j].append(0)

    A = []
    B = []
    for z in range ( int(N) ) :
        a = []
        A.append(a)
        B.append(0)
        for m in range ( int(N) ) :
            A[z].append(0)

    for j in range ( 1 , int(s) + 1 ) :
        A[0][0] = - 3 * alpha / 4 - 3
        A[0][1] = 1
        B[0] = - 3 * alpha / 4 * P[j-1][0] - 2 * P_L
        for z in range ( 1 , int(N-1) ) :
            A[z][z-1] = A[z][z+1]  = 1
            A[z][z] = - alpha - 2
            B[z] = - alpha * P[j-1][z]

        A[int(N)-1][int(N)-2] = 1
        A[int(N)-1][int(N)-1] = - 3 * alpha / 4 - 3
        B[int(N)-1] = - 3 * alpha / 4 * P[j-1][int(N)-1] - 2 * P_R

        #P[j] = list(la.solve(A , B))
        A = np.array(A)
        B = np.array(B)
        P[j] = list(np.linalg.inv(A).dot(B))

    P_list_im[i] = P

for i  in range ( len ( delta_x ) ) :
    for j in range ( int ( np.log10(S[i]) ) ) :
        for q in [ 1 , 3 , 5 , 8 ] :
            z = q * 10**j
            P_list_im[i][z] = [P_L] + P_list_im[i][z] + [P_R]
            plt.figure('implicit - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.title('implicit - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.xlabel('x(m)')
            plt.ylabel('Pressure(pa)')
            plt.plot( X_list[i] , P_list_im[i][z] , label = '∆t = {}'.format( T_list[i][z] ) )
            plt.legend()


# Crank - Nicholson method
P_1_CN = []
P_2_CN = []
P_3_CN = []
P_list_CN = [ P_1_CN , P_2_CN , P_3_CN ]

for i in range ( len ( delta_x ) ) :
    dx = delta_x[i]
    dt = delta_t[i]
    N = L / dx
    P = P_list_CN[i]
    s = S[i]
    beta = alpha / dt * dx**2

    for j in range ( int(s) + 1 ) :
        a = []
        P.append(a)
        for j1 in range ( int(N) ) :
            P[j].append(0)

    A = []
    B = []
    for z in range ( int(N) ) :
        a = []
        A.append(a)
        B.append(0)
        for m in range ( int(N) ) :
            A[z].append(0)

    for j in range ( 1 , int(s) + 1 ) :
        A[0][0] = -3 * beta / 2 - 3
        A[0][1] = 1
        B[0] = ( -3 * beta / 2 + 3 ) * P[j-1][0] - P[j-1][1] - 4 * P_L

        for z in range ( 1 , int(N-1) ) :
            A[z][z-1] = A[z][z+1]  = 1
            A[z][z] = - 2 * beta - 2
            B[z] = -P[j-1][z-1] + ( -2 * beta + 2 ) * P[j-1][z] - P[j-1][z+1]

        A[int(N)-1][int(N)-2] = 1
        A[int(N)-1][int(N)-1] = - 3 * beta / 2 - 3
        B[int(N)-1] = - P[j-1][int(N)-2] + ( -3 * beta / 2 + 3 ) * P[j-1][int(N)-1] - 4 * P_R

        P[j] = list(la.solve(A , B))
    P_list_CN[i] = P

for i  in range ( len ( delta_x ) ) :
    for j in range ( int (np.log10(S[i]) ) ) :
        for q in [ 1 , 3 , 5 , 8 ] :
            z = q * 10**j
            P_list_CN[i][z] = [P_L] + P_list_CN[i][z] + [P_R]
            plt.figure('Crank-Nicholson - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.title('Crank-Nicholson - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.xlabel('x(m)')
            plt.ylabel('Pressure(pa)')
            plt.plot( X_list[i] , P_list_CN[i][z] , label = '∆t = {}'.format( T_list[i][z] ) )
            plt.legend()


## Analytical method
P_1_A = []
P_2_A = []
P_3_A = []
P_list_A = [ P_1_A , P_2_A , P_3_A ]

for i in range ( len ( delta_x ) ) :
    dx = delta_x[i]
    dt = delta_t[i]
    N = L / dx
    P = P_list_A[i]
    s = S[i]
    T = T_list[i]
    X = X_list[i]

    for j in range ( int(s) + 1 ) :
        a = []
        P.append(a)
        for j1 in range ( int(N) ) :
            P[j].append(0)

    for t in range( 1 , int(s) + 1 ) :
        for x in range ( int(N) ) :
            n = 1
            sum1 = 0
            while ( ( 1 / n ) * np.exp( - ( n * np.pi / L )**2 * T[t] / alpha ) * np.sin( n * np.pi * X[x+1] / L ) ) > 0.00001 :
                sum1 += ( 1 / n ) * np.exp( - ( n * np.pi / L )**2 * T[t] / alpha ) * np.sin( n * np.pi * X[x+1] / L )
                n += 1

            P[t][x] = P_L + ( P_R - P_L ) * ( X[x+1] / L + 2 / np.pi * sum1 )
            if P[t][x] < 0 :
                P[t][x] = 0

    P_list_A[i] = P

for i  in range ( len ( delta_x ) ) :
    for j in range ( int (np.log10(S[i]) ) ) :
        for q in [ 1 , 3 , 5 , 8 ] :
            z = q * 10**j
            P_list_A[i][z] = [P_L] + P_list_A[i][z] + [P_R]
            plt.figure('Analytical - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.title('Analytical - ∆x = {}cm'.format( int ( 100 * delta_x[i] ) ) )
            plt.xlabel('x(m)')
            plt.ylabel('Pressure(pa)')
            plt.plot( X_list[i] , P_list_A[i][z] , label = '∆t = {}'.format( T_list[i][z] ) )
            plt.legend()
plt.show()







