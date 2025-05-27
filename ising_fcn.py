import numpy as np
import math as m
import matplotlib.pyplot as plt
import random as r



def matrice0(n):
    A = np.random.choice([-1, 1], size = (n, n))

    return (A)

def matr(A):
    n=len(A)
    B = np.vstack((A[1:,:], A[0,:])) #moves up (lower neighbour)
    C = np.vstack((A[-1,:], A[0:-1,:])) #moves down (upper neighbour)
    D = np.hstack((A[:,1:], np.reshape(A[:,0],(n,1)))) #moves left (neighbour to the right)
    E = np.hstack((np.reshape(A[:,-1],(n,1)), A[:,0:-1])) #moves right (neighbour to the left)
  
    return (B, C, D, E)
  

def calculEandshift(A, T, N, J, H):
    Kb=1
    B, C, D, E =matr(A)
    i = r.randint(0, N-1)
    j = r.randint(0, N-1)
    deltaE= 2*J*A[i][j]*(B[i][j]+C[i][j]+D[i][j]+E[i][j]) + H*A[i][j]
    
    p = m.exp(-deltaE/(Kb*T))
    
    
    w = r.random()
    
    if p>w:
        A[i][j] = -A[i][j]

    return A

#============================================FUNCTION ISING========================================
def ising(A,numiter,N,T,H, graph):
    "the ising function"
    J = 1.
    k = 1.
    #generate empy arrays for storing Momentum M and energy E
    M=np.array([]);
    E=np.array([]);
    # Generate initial spins matrix
    #.....
    #Initialize mean value, in time, of the Momentum and Energy (i.e. "measured quantities")
    if graph == "yes":
        plt.figure(1)
        plt.show()
    #plt.subplots(figsize=(10, 20))
    grid = A
    for i in range(numiter+1): 
        B, C, D, F = matr(grid)
        Mi = np.sum(grid)
        Ei = np.sum(-0.5*J*grid*(B+C+D+F) - H*grid) 
        M = np.append(M, Mi)
        E = np.append(E, Ei)
        grid = calculEandshift(grid, T, N, J, H)
        
        #the array of spins is called grid, filled with +1 and -1, representing spins up and down
        #print every 1000 iteration, in this example
        if (i%1000 == 0 and graph == "yes"):
                plt.subplot(4,1,1)
                plt.xlabel("i= %7d, T = %0.2f, M = %0.2f, E = %0.2f" %(i, T, M[i], E[i]))
                plt.imshow((grid+1)*128,cmap='bone',interpolation="nearest");
                
                plt.subplot(4,1,3)
                tgrid = range(1,i+1);
                plt.plot(tgrid, M[:-1],'bo')
                plt.xlabel("iteration")
                plt.ylabel("Inst Mom")
                #plt.show()
                
                #plt.figure(3)
                plt.subplot(4,1,4)
                tgrid = range(1,i+1);
                plt.plot(tgrid, E[:-1],'go')
                plt.xlabel("iteration")
                plt.ylabel("Inst En")
                plt.show() 
#                plt.pause(2)
                plt.pause(0.00000001)
    #mean values in time
    Mmean= Mi/(N**2)
    Emean= Ei/(N**2)
    Mmean2=np.sum(grid**2)/(N**2)
    Emean2=np.sum((-0.5*J*grid*(B+C+D+F) - H*grid)**2) /(N**2)        
    return Mmean,Emean,Mmean2,Emean2




def ising_temp(iterations, N):
    
    
    
    MmeanT=np.array([]);
    EmeanT=np.array([]);
    #MmeanT2=np.array([]);
    #EmeanT2=np.array([]);
    chi=np.array([]);
    cv=np.array([]);
    
    
    grid = np.ones((N, N))
    Ti = np.arange(0.1, 4.1, 0.1)
    H = 0
    
    for T in Ti:
        Mmeani,Emeani, Mmeani2, Emeani2 = ising(grid, iterations, N, T, H, graph = "no")
        MmeanT = np.append(MmeanT, Mmeani)
        EmeanT = np.append(EmeanT, Emeani)
        #MmeanT2 = np.append(MmeanT2, Mmeani2)
        #EmeanT2 = np.append(EmeanT2, Emeani2)
        
        chi_i = 1/T*(Mmeani2 - Mmeani**2)
        cv_i = 1/T**2*(Emeani2 - Emeani**2)
        
        chi = np.append(chi, chi_i)
        cv = np.append(cv, cv_i)
        
        percentage = int(T/Ti[-1]*100)
        
        print("%d%%" % percentage)
    
    plt.figure(1)
    
    #plot chi or EmeanT
    plt.subplot(231)
    #plt.plot(Ti, EmeanT, 'ro');
    #plt.ylabel('energy per site');
    plt.plot(Ti, chi, 'ro');
    plt.ylabel('susceptibility');
    plt.xlabel('Temperature');
    
    #plot cv or MmeanT
    plt.subplot(233)
    #plt.plot(Ti, MmeanT, 'bo');
    #plt.ylabel('Magnetization per site');
    plt.plot(Ti, cv, 'bo');
    plt.ylabel('cv');
    plt.xlabel('Temperature');
    plt.show()
    
    #plot chi or EmeanT
    plt.subplot(234)
    plt.plot(Ti, EmeanT, 'ro');
    plt.ylabel('energy per site');
    plt.xlabel('Temperature');
    
    #plot cv or MmeanT
    plt.subplot(236)
    plt.plot(Ti, MmeanT, 'bo');
    plt.ylabel('Magnetization per site');
    plt.xlabel('Temperature');
    plt.show()
    
    return




def ising_field(iterations, N):
    
    chi1=np.array([]);
    cv1=np.array([]);
    chi2=np.array([]);
    cv2=np.array([]);
    
    grid=matrice0(N)
    Ti = (4.0, 8.0)
    Hi = (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
    counter = 0
    
    for T in Ti:
        for H in Hi:
            counter = counter + 1
            Mmeani,Emeani, Mmeani2, Emeani2 = ising(grid, iterations, N, T, H, graph = "no")
            
            chi_i = 1/T*(Mmeani2 - Mmeani**2)
            cv_i = 1/T**2*(Emeani2 - Emeani**2)
            
            if T == Ti[1]:
                chi1 = np.append(chi1, chi_i)
                cv1 = np.append(cv1, cv_i)
            else:
                chi2 = np.append(chi2, chi_i)
                cv2 = np.append(cv2, cv_i)
                
            percentage = int(counter/(len(Hi) * len(Ti)) * 100)
            print("%d%%" % percentage)
            
            
    plt.figure(1)
    plt.subplot(231)
    #plt.plot(Ti, EmeanT, 'ro');
    #plt.ylabel('energy per site');
    plt.plot(Hi, chi1, 'ro');
    plt.ylabel('susceptibility (T=4.0)');
    plt.xlabel('Field H');
    
    #plot cv or MmeanT
    plt.subplot(233)
    #plt.plot(Ti, MmeanT, 'bo');
    #plt.ylabel('Magnetization per site');
    plt.plot(Hi, cv1, 'bo');
    plt.ylabel('cv (T=4.0)');
    plt.xlabel('Field H');
    
    plt.subplot(234)
    #plt.plot(Ti, EmeanT, 'ro');
    #plt.ylabel('energy per site');
    plt.plot(Hi, chi2, 'ro');
    plt.ylabel('susceptibility (T=8.0)');
    plt.xlabel('Field H');
    
    #plot cv or MmeanT
    plt.subplot(236)
    #plt.plot(Ti, MmeanT, 'bo');
    #plt.ylabel('Magnetization per site');
    plt.plot(Hi, cv2, 'bo');
    plt.ylabel('cv (T=8.0)');
    plt.xlabel('Field H');
    plt.show()
    
    return