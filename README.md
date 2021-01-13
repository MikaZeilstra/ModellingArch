### Modelling archearhodopsion
This is a project to model the photocycle of NovArch by fitting data to potential models. Therefore this python module was created. It is meant to be fairly easy to use and applicable in any situation where there are transitions which are linearly dependant on an outside stimulus.

### Dependencies
The tool is designed using Python 3.9

The required packages are :
- matplotlib
- sympy
- scipy

### Usage on its own
Firstly download the package using pip:
```
pip install ModellingArch
```
Now you can use the tool using the following command:
```
python -m ModellingArch -h
```

An example using some of the options is:
```
python -m ModellingArch -Matrix "[[-k0,0,k2],[k0,-k1,0],[0,k1,-k2]]" -Ft "[(1,k1)]" -Lt "[k0]" -LightIntensities "[1,0.5]" -Data "[[10],[5]]" -FixedT "[(k1,0.1)]"
```

### Usage as a module
This package can ofcourse also be used as a module in your own python program. It provides a single class Model which contains everything you need to fit the data to a model.
