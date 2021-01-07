import sympy as sp
import scipy.optimize as opt
import scipy.integrate as inte
import matplotlib.pyplot as plt
import numpy as np
sp.init_printing(use_unicode=False, wrap_line=False)

Manualinput = True

#The index of the state from which the fluorecent transition comes and the index of the fluorecent transitions as tuples.
FluorecentStatesAndTransitions = [(2,3),(2,4)]


#Maximum number of tries to solve for the transitions until error, Interval in which the intialguess will be randomly chosen and which rates are considerd to be 0
MaxTries = 50
IntialGuessInterval = (0,100)
ZeroEpsilon = 1e-4

#Setup configuration for numerical integration Firststepsize is for accuracy at the start and maxstepsize influences the resolution
SolutionInterval = (0,10)
FirstStepSize = 1e-4
MaxStepSize = 1e-1

#Setup Use of Predetermined Rates
UsePredeterminedRates = False
PreRates = [10,10,10,1]

#Fixed transition rates can be inserted here as transition-value tuples THEY NEED TO BE SORTED IN ASCENDING ORDER OF TRANSITION INDEX
FixedK = []


#Initiate the T symbol
t = sp.symbols("t", positive=True,real=True)

if Manualinput :
    # Fill in the number of transitions and states
    NumberOfTransitions = 7
    NumberOfStates = 4

    #Fill in the light intensities, and the affected transitions
    IntensityFractions = [0.95,0.75,0.5,0.25]
    LightDependantTransitions = [0,2,5]

    #Fill in the K_totals corresponding to the given light intensities.
    KTotals = [[-0.6958199483,-400],[-0.6062190023,-300],[-0.530258567,-200],[-0.447416794,-100]]

    #Initialize the symbols for the transitions, states and the time
    k = []
    for Kn in range(NumberOfTransitions):
        k.append(sp.symbols("k" + str(Kn), positive=True,real=True))



    #Setup matrix using the just created symbols
    M = sp.Matrix([[-k[0] ,   k[2],0,k[6]],
                    [k[0]   ,-k[2]-k[1],k[3],0],
                    [0,k[1],-k[3]-k[4],k[5]],
                    [0,0,k[4],-k[5]-k[6]]])




#A lot of The manual features are not in the CLI version
else:
    raise NotImplementedError()
    # NumberOfStates = int(input("Enter the number of states:"))
    #
    # # Initialize matrix
    # InputMatrix = []
    # print("Enter the entries of the transition matrix rowwise:")
    #
    # # For user input
    # for i in range(NumberOfStates):  # A for loop for row entries
    #     a = []
    #     for j in range(NumberOfStates):  # A for loop for column entries
    #         a.append(input())
    #     InputMatrix.append(a)
    #
    # # For printing the matrix
    # #print(InputMatrix)
    #
    # #Convert matrix to sympy
    # M = sp.Matrix(sp.sympify(InputMatrix))
    # k = list(M.free_symbols)
    # print(k)
    #
    # if(len(k) != NumberOfTransitions):
    #     raise RuntimeError("Too many symbols in matrix for number of transitions")
    # #M = sp.Matrix(M)
    #
    # # Fill in the known lambdas and number of transitions
    # NumberOfLambdas = int(input("Enter the number of lambdas:"))
    #
    # print("Enter the lamdas:")
    # KnownLambda = []
    # for j in range(NumberOfLambdas):  # A for loop for column entries
    #     KnownLambda.append(int(input()))
    # print(KnownLambda)



#Setup Symbols for the states
s = []
for Sn in range(NumberOfStates):
    s.append(sp.symbols("s" + str(Sn), positive=True, real=True))

#Skip calculating rates if we were given rates
if (not UsePredeterminedRates):
    #Find the charactertic polynomial
    cPoly = M.charpoly().expr

    #Get the lambda symbol
    lam = cPoly.free_symbols.difference(sp.FiniteSet(*k)).pop()

    #Substitute the lambdas for the known values and scale affected transition rates to the given value
    filledLambdas = []
    for i in range(len(IntensityFractions)):
        dependentk = list(map(lambda x: k[x], LightDependantTransitions))
        ExtraCPoly = cPoly.subs(zip(dependentk, map(lambda x: IntensityFractions[i] * x, dependentk)))
        filledLambdas += list(map(lambda x: ExtraCPoly.subs(lam, x), KTotals[i]))

    #Store the actual symbol instead of reference in Fixed ks
    FixedK = list(map(lambda x: (k[x[0]], x[1]), FixedK))
    for i in range(len(filledLambdas)):
        filledLambdas[i] = filledLambdas[i].subs(FixedK)


    #Use set subtraction to get the symbols for which need to be solved
    UnkownK = [k for k in k if k not in map(lambda k : k[0], FixedK)]
    NumberOfUnkownTransitions = len(UnkownK)


    #Export system to list of functions
    def func(fun):
        return lambda x : fun(*x)
    NonlinearSystem = []
    for exp in filledLambdas:
        NonlinearSystem.append(func(sp.lambdify(UnkownK,exp)))


    #Setup the variables needed to loop
    Tries = 0
    Rates = {"success" : False,"active_mask":[-1]}

    #Try and find a solution with new random starting points 10 times. a solution is reject if the bounds are active (active_mask) this means the solution is either on or lower than th given lower bound
    while (((not(all([x == 0 for x in Rates['active_mask']]))) or (not Rates["success"])) and (Tries < MaxTries)):
        Rates = opt.least_squares(lambda num : list(map(lambda x: x(num), NonlinearSystem)), np.random.rand((NumberOfUnkownTransitions))*(IntialGuessInterval[1]-IntialGuessInterval[0])+IntialGuessInterval[0],bounds=(ZeroEpsilon,np.inf))
        Tries += 1

    #Print the solution if found otherwise error
    if(Tries < MaxTries):
        print("Found transition rates in " + str(Tries) + " tries.")
        #The k values are the value of x in the Rates object. These values are the free transition rates in order with the known transition rates removed
        print(Rates)
    else:
        #If we havenot find a solution throw an error and show the last Rates Object.
        print(Rates)
        raise RuntimeError("Could not solve for the transition rates in " + str(MaxTries)+ " tries with given constraints Not all in bounds was " + str((not(all([x == 0 for x in Rates['active_mask']])))))
else:
    Rates = {"x" : PreRates}

#Setup list of the transition rates
TransitionRates = list(Rates["x"])

#Add the Fixed Transtion Rates back in
for x in FixedK:
    TransitionRates.insert(k.index(x[0]),x[1])

print("Found Transition Rates " + str(TransitionRates))


#Multiply the light dependend transition with a pattern
#TransitionRates[0] *= sp.Abs(sp.cos(t))

#Setup system of ODEs using matrix vector multiplication
ODEs = M.subs(zip(k,TransitionRates)) * sp.Matrix(s)

#Convert the system to a function for use by scypy
ODESystem = lambda ti,y : list(map(lambda x : x[0], sp.lambdify([t] + s,ODEs)(ti,*y)))


#solve the system
solution = inte.solve_ivp(ODESystem,SolutionInterval,[1] + [0]* (NumberOfStates-1),first_step= FirstStepSize,max_step = MaxStepSize,method="Radau")

#Plot the results
for i in range(len(solution['y'])):
    plt.plot(solution["t"], solution["y"][i], label="State " + str(i))

#Plot the fluorecence as the population which transitions from S1 with Rate k1 for all the given fluorecent transitions
for fluorecence in FluorecentStatesAndTransitions:
    plt.plot(solution["t"],solution["y"][fluorecence[0]]*TransitionRates[fluorecence[1]],label="Fluorecense Transition " + str(fluorecence[1]))
plt.legend()
plt.show()







