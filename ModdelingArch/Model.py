import sympy as sp
import scipy.optimize as opt
import scipy.integrate as inte
import numpy as np
import matplotlib.pyplot as plt


class Model:
    M = sp.Matrix
    LightDependantTransitions = list[sp.Symbol]
    FluorescentTransitions = list[tuple[sp.Symbol,int]]

    NumberOfStates = int
    States = list[sp.Symbol]

    NumberOfTransitions = int
    Transitions = list[sp.Symbol]

    TransitionRates = dict[sp.Symbol : float]

    SolutionTimes = list[float]
    Solution = list[list[float]]

    '''
    Constructor for the model class
    :param Matrix : Matrix describing the state diagram of the continuous markov chain model containing free sympy variables for each transition
    :type Matrix : sympy.Matrix
    :param LightDependentTransitions : List of symbols coresponding to the  dependent transitions of the model
    :type LightDependentTransitions : list[sympy.Symbol]
    :param FluorescentTransitions : List of tuples with the first value being a integer representing the state number (0- indexed) which the transition originates from and the second value being the symbol it orginates from
    :type FluorescentTransitions : list[tuple[sp.Symbol,int]]
    '''
    #Override
    def __init__(self, Matrix : sp.Matrix, LightDependentTransitions :list[sp.Symbol], FluorescentTransitions : list[tuple[sp.Symbol,int]]):
        self.M = Matrix
        self.LightDependantTransitions = LightDependentTransitions
        self.FluorescentTransitions = FluorescentTransitions

        self.NumberOfStates = self.M.shape[0]
        self.States = list(map(lambda x :  sp.symbols("s" + str(x), positive=True, real=True), range(self.NumberOfStates)))

        self.Transitions = list(self.M.free_symbols)
        self.NumberOfTransitions = len(self.Transitions)

    '''
    A function to find the transition rates of a model given data for different intensities. Stores the found values in the model and returns them.
    :param data : Dictionary containing measured eigenvalues (negative apparent rates) in a list as values for intensities as keys
    :type data : dict[float, list[float]]
    :param FixedTransitions : Dictionary containing symbols of transitions as keys and their values, defaults to {}
    :type FixedTransitions : dict[sympy.Symbol, float]
    :param InitialGuessInterval : The interval from which the random guesses for the transition rates will be taken , defaults to (0,10)
    :type InitialGuessInterval : tuple[float,float]
    :param LowerBoundTransitionRates : The lowest transition rate we allow, defaults to 1e-10
    :type LowerBoundTransitionRates : float
    :param MaxTries: maximum amount of times a new random initial guess is tried, defaults to 50
    :raise RuntimeError : if there is no valid solution for the given data
    :return: A dictionary containing the transitions as keys and their rates as values
    :rtype: dict[sp.symbol,float]
    '''
    def calculate_transitions(self, data : dict[float, list[float]], FixedTransitions=None, InitialGuessInterval= (0,10),LowerBoundTransitionRates = 1e-10,MaxTries=50 ) -> dict[sp.Symbol,float]:
        if FixedTransitions is None:
            FixedTransitions = {}

        # Find the charactertic polynomial
        cPoly = self.M.charpoly().expr

        # Get the lambda symbol
        lam = cPoly.free_symbols.difference(sp.FiniteSet(*self.Transitions)).pop()

        # Substitute the lambdas for the known values and scale affected transition rates to the given value
        filledLambdas = []
        for Intensity in data.keys():
            ExtraCPoly = cPoly.subs(zip(self.LightDependantTransitions, map(lambda x: Intensity * x, self.LightDependantTransitions)))
            filledLambdas += list(map(lambda x: ExtraCPoly.subs(lam, x), data[Intensity]))

        # Store the actual symbol instead of reference in Fixed ks
        for i in range(len(filledLambdas)):
            filledLambdas[i] = filledLambdas[i].subs(FixedTransitions.items())
        
        
        # Use set subtraction to get the symbols for which need to be solved
        UnkownK = [k for k in self.Transitions if k not in FixedTransitions.keys()]
        NumberOfUnkownTransitions = len(UnkownK)

        # Export system to list of functions
        def func(fun):
            return lambda x: fun(*x)

        NonlinearSystem = []
        for exp in filledLambdas:
            NonlinearSystem.append(func(sp.lambdify(UnkownK, exp)))

        # Setup the variables needed to loop
        Tries = 0
        Rates = {"success": False, "active_mask": [-1]}

        # Try and find a solution with new random starting points 10 times. a solution is reject if the bounds are active (active_mask) this means the solution is either on or lower than th given lower bound
        while (((not (all([x == 0 for x in Rates['active_mask']]))) or (not Rates["success"])) and (Tries < MaxTries)):
            Rates = opt.least_squares(fun=lambda num: list(map(lambda x: x(num), NonlinearSystem)),
                                      x0=np.random.rand((NumberOfUnkownTransitions)) * (InitialGuessInterval[1] - InitialGuessInterval[0]) +InitialGuessInterval[0],
                                      bounds=(LowerBoundTransitionRates, np.inf))
            Tries += 1

        # Print the solution if found otherwise error
        if (Tries >= MaxTries):
            # If we havenot find a solution throw an error and show the last Rates Object.
            print(Rates)
            raise RuntimeError("Could not solve for the transition rates in " + str(
                MaxTries) + " tries with given constraints Not all in bounds was " + str(
                (not (all([x == 0 for x in Rates['active_mask']])))))

        print("Found transition rates in " + str(Tries) + " tries.")
        # The k values are the value of x in the Rates object. These values are the free transition rates in order with the known transition rates removed
        print(Rates)

        # Setup list of the transition rates
        self.TransitionRates = dict(zip(UnkownK, Rates['x'])) | FixedTransitions

        return self.TransitionRates

    '''
    Solves the differential system of equations for the population of each state as a function of time and shows the results
    :param SolutionInterval
    '''
    def find_population(self,SolutionInterval =(0,10),FirstStepSize=1e-4, MaxStepSize =1e-1):
        ODEs = self.M.subs(self.TransitionRates.items()) * sp.Matrix(self.States)

        # Convert the system to a function for use by scypy
        ODESystem = lambda ti, y: list(map(lambda x: x[0], sp.lambdify([sp.Symbol("t")] + self.States, ODEs)(ti, *y)))

        # solve the system
        solution = inte.solve_ivp(ODESystem, SolutionInterval, [1, 0, 0], first_step=FirstStepSize,
                                  max_step=MaxStepSize, method="Radau")



        # Plot the results
        for i in range(len(solution['y'])):
            plt.plot(solution["t"], solution["y"][i], label="State " + self.States[i].name[1:])

        # Plot the fluorecence as the population which transitions from S1 with Rate k1 for all the given fluorecent transitions
        for fluorecence in self.FluorescentTransitions:
            plt.plot(solution["t"], solution["y"][fluorecence[0]] * self.TransitionRates[fluorecence[1]],
                     label="Fluorecense Transition " + str(fluorecence[1].name))
        plt.legend()
        plt.show()


