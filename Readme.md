### Moddeling archearhodopsion
This is a project to model the photocycle of NovArch by fitting data to potential models. The program in this repository serves as a tool to automate this process.

### Requirements
The tool is designed using Python 3.9

The required packages are :
- matplotlib
- sympy
- scipy

### Usage
First The number of states and transitions is inserted into their respective variables. After that the matrix representing the model needs to be inserted at the M variable. The fluorescent states and transitions can be inserted at the FluorecentStatesAndTransitions variable. The light dependant transitions can be defined at its respective variable as well.

When the model is inserted using the above method, The measurements need to be added.
First a list of the intensities used to measure the eigenvalues needs to be given in the IntensityFractions variable. Lastly the actual (negative) eigenvalues need to be inserted in a 2 dimensional list containing the measurements for each intensity in the same order.

At the start of the file several other variables are also provided to configure the tool more precisely.

After this the program should run print the found true transition rates in the console and show a graph of the state populations and fluorescence.