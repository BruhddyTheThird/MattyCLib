import sympy as sp
import os.path
import sys

# go up one directory level from this file's directory:
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# prepend parent directory to the system path:
sys.path.insert(0, path)


import FunctionLib.MTHFunctions as mtf
import FunctionLib.ParametricGraphing as pmg
import matplotlib.pyplot as plt
import numpy as np
x, y, z, u, v, w, s, t= sp.symbols('x,y,z,u,v,w,s,t')
a,b = sp.symbols('a,b',)
r, theta = sp.symbols('r, theta',)
p, phi = sp.symbols('p, phi')
pi = sp.pi
sin, cos, tan, sqrt, e= sp.sin, sp.cos, sp.tan, sp.sqrt, sp.exp

Field1 = sp.Array([
    7*x+y**2,
    z**2-2*y,
    z+x**2
])
Surface1 = sp.Array([
    sp.cos(t),
    0,
    -sp.sin(t)
])
Surface2 = sp.Array([
    s*sp.cos(t),
    4,
    s*sp.sin(t)
])
Expr1 = 3*(z*x**2+z*y**2)

Point1 = [7,0,0]


parameterField1 = [
    y,
    z,
]


t = sp.Symbol('t')

parameterField = [
    r*sp.cos(theta),
    r*sp.sin(theta)
]


varList1 = [x,y,z]
varList2 = [x,y]


varTupleList1 = [(x,y)]

parVarList1 = [t,s]

parBoundsTupleList1 = [(0,1),(0,1)]
parBoundsTupleList2 = [(0,2*sp.pi)]

div1 = mtf.Div(Field1)
ediv1 = mtf.Div(Field1,Point1)
print('Divergence expression is:\n',div1,'\nDivergence evaluation is:\n',ediv1)
curl1 = mtf.Curl(Field1,VarList=varList1)
ecurl1 = mtf.Curl(Field1,Point1,varList1)
print('Curl expression is:\n',curl1,'\nCurl evaluation is:\n',ecurl1)

#try:
#    potFun = findPotntialFnction(Field1,VarList=varList1)
#    print('A potential function is: \n',potFun)
#except:
#    print('Field:\n',Field1,'\nIs not conservative. Therefore a potential function does not exist for the field.')
Oint = mtf.CountourInt(Field1,varList1,VarBounds=parBoundsTupleList1,Eval=True,Hard=False,ParVarBounds=parBoundsTupleList2,Reversed=False)
print(Oint)
#intExp1 = (mtf.ParametrizeExpr(Expr1,Surface1,VarList=varList1))
#print(mtf.SurfArea(Surface1,varList2,ParVarList=parVarList1,ParaField=parameterField1,Eval=True,Easy=False,ParVarBounds=parBoundsTupleList1,Left=False,IntExpr=intExp1))
#print(mtf.Magnitude(mtf.CrossProduct([-v*sp.sin(u),v*sp.cos(u),1],[sp.cos(u),sp.sin(u),0])))
#NormalPre = []
#for j in range(2):
#    NormalPre.append([sp.trigsimp(sp.diff([v*sp.cos(u),v*sp.sin(u),u][i],[u,v][j])) for i in range(3)])
#print(NormalPre)
#pmg.NumpyParaGraph(Field1,varTupleList1,parBoundsTupleList1)
#plt.show()
#Int1 = mtf.SurfFlux(Surface1,Field1,varList1,parVarList1,parBoundsTupleList1,Eval=True,Left=False,Upward=True)
#print(Int1)
