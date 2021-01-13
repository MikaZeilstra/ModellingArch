import argparse
import ast
import sympy as sp
import numpy as np
from ModellingArch.Model import Model

parser = argparse.ArgumentParser(
                                 description="A pre-made tool to solve the problem of fitting data to a coninuus markov procces with a given model.\n\n"+
                                            "Pick either an input file or insert the model using the arguments, for larger models it is recommended to mannually make an input file",
                                 allow_abbrev=False)
#Model
parser.add_argument('-Matrix',dest="M", required=False, help ="The square matrix describing the markov chain state diagram where each row will give the diffential equation for the population of that state")
parser.add_argument('-Ft','-FluorescentTransitions',dest="Ft", required=False, help="A list containing the names of each fluorescent transition")
parser.add_argument('-Lt','-LightDependentTransitions', dest="Lt",required=False, help="A list containing the names of each light dependent transition")

#Data
parser.add_argument('-Li','-LightIntensities',dest="Li", required=False, help="The light intensities the data was gathered at")
parser.add_argument('-Data', required=False, dest="Data",help= "A 2-dimensional list containing a set of measured apparent rates for each light intensity")

#FixedTransitions
parser.add_argument('-FixedT -FixedTransitions', required=False, dest="FixedT", help="A list of ordered pairs containing the name of the transition and its value for each fixed transition", default="[]")

#Transition finding settings
parser.add_argument("-Ig","-InitialGuessInterval", dest="Ig", required=False, help="A tuple containing the interval from which the initial guesses will be taken",default="(0,10)")
parser.add_argument("-L","-TransitionsLowerBound", dest="Lb", required=False, help="A float which determines the lowest accepted transition rate" , default="10e-10")
parser.add_argument("-Mt","-MaxTries", dest="Mt", required=False, help="An integer detemining how many times we want to try and find the transition rates if it has failed with current intial values", default="50")

#Solution finding settings
parser.add_argument("-Si","-SolutionInterval", dest="Si", required=False, help="A tuple containing the interval for which the populations will be calculated, defaults to (0,10)",default="(0,10)")
parser.add_argument("-Ip","-InitialPopulation", dest="Ip", required=False, help="A list containing the initial populations of the system", default="None")
parser.add_argument("-Fs","-FirstStepSize", dest="Fs", required=False, help="The size of the first step taken while numerically solving the system the next step will be dynamically determined", default="1e-4")
parser.add_argument("-Ms","-MaxStepSize", dest="Ms", required=False, help="The maximum stepsize that will be taken to calculate the populations determines the resolution of the solution", default="1e-1")
args = parser.parse_args()

#Check if the model is given using the command line or a file
if (not(args.M is None or args.LightIntensities is None or args.Data is None or args.Lt is None or args.Ft is None)):
    M = sp.Matrix(sp.sympify(args.M))

    k = sorted(list(M.free_symbols), key=lambda x : x.name)

    LightIntensities = ast.literal_eval(args.Li)
    LightData = ast.literal_eval(args.Data)

    Ldt = sp.sympify(args.Lt)
    Ft = sp.sympify(args.Ft)

    Data = None
    if len(LightData) != len(LightIntensities):
        raise Exception("Length of LightIntensities and Data must be the same")
    else:
        Data = dict(zip(LightIntensities,np.array(LightData) * -1))
else:
    raise IOError("Please input all of MATRIX ,FT ,LT , LI and DATA, Or give INPUTFILE")

#Import the fixed transitions
FixedT = dict(sp.sympify(args.FixedT))

#Import the options for the calculating of transitions
InitialGuess = ast.literal_eval(args.Ig)
LowerBound = ast.literal_eval(args.Lb)
MaxTries = ast.literal_eval(args.Mt)

#Import the options for calculating the populations
SolutionInterval = ast.literal_eval(args.Si)
InitialValue = ast.literal_eval(args.Ip)
FirstStep = ast.literal_eval(args.Fs)
MaxStep = ast.literal_eval(args.Ms)

#Run the model
model = Model(M,Ldt,Ft)
model.calculate_transitions(data=Data,FixedTransitions=FixedT, InitialGuessInterval=InitialGuess,LowerBoundTransitionRates=LowerBound, MaxTries=MaxTries)
model.find_population(SolutionInterval=SolutionInterval,InitialCondition=InitialValue,FirstStepSize=FirstStep,MaxStepSize=MaxStep)