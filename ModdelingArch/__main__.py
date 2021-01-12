import argparse
import ast
import sympy as sp
from Model import Model

parser = argparse.ArgumentParser(description="A pre-made tool to solve the problem of fitting data to a coninuus markov procces with a given model.\n\n"+
                                            "Pick either an input file or insert the model using the arguments, for larger models it is recommended to mannually make an input file")
parser.add_argument('-Matrix',dest="M", required=False, help ="The square matrix describing the markov chain state diagram where each row will give the diffential equation for the population of that state")
parser.add_argument('-Ft',dest="Ft", required=False, help="A list containing the names of each fluorescent transition")
parser.add_argument('-Lt', dest="Lt",required=False, help="A list containing the names of each light dependent transition")
parser.add_argument('-LightIntensities',dest="LightIntensities", required=False, help="The light intensities the data was gathered at")
parser.add_argument('-Data', required=False, dest="Data",help= "A 2-dimensional list containing a set of measured apparent rates for each light intensity")

args = parser.parse_args()

M = sp.Matrix(sp.sympify(args.M))

k = sorted(list(M.free_symbols), key=lambda x : x.name)

km = sp.symbols("k0 k1 k2 k3")

LightIntensities = ast.literal_eval(args.LightIntensities)
LightData = ast.literal_eval(args.Data)

Ldt = sp.sympify(args.Lt)
Ft = sp.sympify(args.Ft)

Data = None
if len(LightData) != len(LightIntensities):
    raise Exception("Length of LightIntensities and Data must be the same")
else:
    Data = dict(zip(LightIntensities,LightData))

model = Model(M,Ldt,Ft)
model.calculate_transitions(data=Data)
model.find_population()