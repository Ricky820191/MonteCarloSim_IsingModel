import ising_fcn as t
import numpy as np
import matplotlib.pyplot as plt


choice = 0
while choice != "4":
    
    print(" ")
    print(" ")
    print("""Welcome, this program uses the ising model to compute several variables
    in various conditions:
        1: Evolution of the system over time starting from a random state
        2: Evolution of the susceptibility and the specific heat capacity 
           over temperature
        3: Evolution of the same variables of 2 but over the applied field 
           for 2 different temperatures
        4: quit""")
    choice = input(">> ")
    
    if choice == "1":
        
        print("How many spins do you want on one side of the square?")
        N = int(input("     >> "))
        
        print("What should the temperature be? (max value should not be above 10, always positive)")
        T = float(input("     >> "))
        
        print("What should the external field be? (can be positive or negative)")
        H = float(input("     >> "))
        
        print("How many iterations?")
        iterations = int(input("     >> "))
        
        grid=t.matrice0(N)
        
            
        Mmean,Emean,Mmean2,Emean2 = t.ising(grid, iterations, N, T, H, graph = "yes")
        print(" ")
        print(" ")
        print("The value of Mmean is:", Mmean)
        print("The value of Emean is:", Emean)
        print("The value of Mmean2 is:", Mmean2)
        print("The value of Emean2 is:", Emean2)
        print(" ")
        print(" ")
    
    #############################################################################
    elif choice == "2":
        print(" ")
        print('''Here, we are going to do a study of how 
        magnetization and energy per site vary with 
        temperature, hence we will vary the temperature
        in an interval of 0.1 and 4 with a step of 0.1
        it will take a while so be patient, at the end
        the graphs will appear. In this simulation
        we will start with all the spins ALLIGNED and 
        zero external field, sit back and relax until 
        the results appear :)''')  
        print(" ")
        print(" ")
        print("How many spins do you want on one side of the square?")
        N = int(input("     >> "))
        
        print("How many iterations? Note that high values will significantly increase simulation time")
        iterations = int(input("     >> "))
        
        t.ising_temp(iterations, N)
        
    ###############################################################################
    elif choice == "3":
        
        print(" ")
        print(" ")
        print("""Next, we will be changing the applied field from 0.1 to 0.4 
        with a step of 0.1, this will be done twice for different temperatures:
        4.0 and 8.0. Note that we will be assuming that the system is initially
        in a rdndom state. Same as before, this might take some time.""")
        print(" ")
        print(" ")
        
        
        print("How many spins do you want on one side of the square?")
        N = int(input("     >> "))
        
        print("How many iterations? Note that high values will significantly increase simulation time")
        iterations = int(input("     >> "))
        
        t.ising_field(iterations, N)
        
    ###############################################################################
    
    elif choice == "4":
        print("Bye!")
    
    else:
        print("please input either 1, 2, 3 or 4")
        print(" ")
        print(" ")




