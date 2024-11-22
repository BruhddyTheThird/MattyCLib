import sympy as sp
from sympy import diff
from sys import argv
import warnings
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
def InteriorProduct(Array1, Array2):
    """
    Computes the interior (dot) product of two arrays (vectors).

    Parameters:
    Array1 (list or sympy.Array): The first vector.
    Array2 (list or sympy.Array): The second vector, which must have the same length as Array1.

    Returns:
    sympy expression: The sum of the pairwise products of corresponding elements from Array1 and Array2.

    Raises:
    ValueError: If the lengths of Array1 and Array2 are not the same.

    Example:
    >>> Array1 = [x, y, z]
    >>> Array2 = [u, v, w]
    >>> InteriorProduct(Array1, Array2)
    x*u + y*v + z*w
    """
    if len(Array1) != len(Array2):
        raise ValueError('The arrays do not have the same length!')
    productExpr = sum([Array1[i]*Array2[i] for i in range(len(Array1))])
    return productExpr
def CrossProduct(Array1, Array2):
    """
    Computes the cross product of two 2D or 3D vectors.

    Parameters:
    Array1 (list or sympy.Array): The first vector.
    Array2 (list or sympy.Array): The second vector, which must have the same length as Array1 and must be of length 2 or 3.

    Returns:
    list: The cross product of the two vectors as a list of symbolic expressions.

    Raises:
    ValueError: If the vectors are not of length 2 or 3.

    Example:
    >>> Array1 = [x, y, z]
    >>> Array2 = [u, v, w]
    >>> CrossProduct(Array1, Array2)
    [y*w - z*v, z*u - x*w, x*v - y*u]
    """
    if len(Array1) != len(Array2):
        raise ValueError('Arrays must have the same length.')

    if len(Array1) == 3:  # 3D
        return [
            sp.trigsimp(Array1[2]*Array2[1] - Array1[1]*Array2[2]),
            sp.trigsimp(Array1[0]*Array2[2] - Array1[2]*Array2[0]),
            sp.trigsimp(Array1[1]*Array2[0] - Array1[0]*Array2[1])
        ]
    elif len(Array1) == 2:  # 2D
        return [0, 0, sp.trigsimp(Array1[0]*Array2[1] - Array1[1]*Array2[0])]
    else:
        raise ValueError("Only 2D or 3D vectors are supported.")
def Magnitude(ArrayLike):
    """
    Computes the magnitude (norm) of a vector.

    Parameters:
    ArrayLike (list or sympy.Array): The vector whose magnitude is to be computed.

    Returns:
    sympy expression: The magnitude (or Euclidean norm) of the vector.

    Example:
    >>> Magnitude([x, y, z])
    sqrt(x**2 + y**2 + z**2)
    """
    ReturnExp = sp.sqrt(sp.simplify(sum([i**2 for i in ArrayLike])))
    return ReturnExp
def Jacobian(MappedField,VarList):
    """
    Computes the determinant of the Jacobian matrix for a given transformation.

    The Jacobian matrix represents the rate of change of each component of the mapped field
    with respect to the transformation variables.

    Parameters:
    MappedField (list of sympy expressions): The functions representing the mapped field (e.g., x, y, z as functions of u, v, w).
    VarList (list of sympy symbols): The variables with respect to which the mapping is defined.

    Returns:
    sympy expression: The determinant of the Jacobian matrix.

    Raises:
    ValueError: If MappedField and VarList do not have the same length.

    Example:
    >>> MappedField = [x**2 + y, y**2 + z]
    >>> VarList = [u, v]
    >>> Jacobian(MappedField, VarList)
    Determinant of the Jacobian matrix
    """
    if len(VarList) != len(MappedField):
        raise ValueError('MappedField and VarList must have the same length.')
    
    # Compute the Jacobian matrix
    jacobian_matrix = sp.Matrix([[diff(f, v) for v in VarList] for f in MappedField])
    return sp.trigsimp(jacobian_matrix.det())

def Div(Field,Point=None,VarList=None):
    """
    Computes the divergence of a vector field.

    The divergence of a vector field F = <P, Q, R> is the sum of the partial derivatives of
    its components with respect to their respective variables.

    Parameters:
    Field (list of sympy expressions): The components of the vector field.
    Point (list of sympy expressions, optional): A specific point at which to evaluate the divergence.
    VarList (list of sympy symbols, optional): The variables with respect to which to compute the divergence. Defaults to [x, y, z].

    Returns:
    sympy expression: The divergence of the vector field, or its value at the specified point.

    Example:
    >>> Field = [x**2, y**2, z**2]
    >>> Div(Field)
    2*x + 2*y + 2*z
    """
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    divExpr = sum([diff(f,s) for f,s in zip(Field,VarList)])
    if Point == None:
        return divExpr
    else:
        for p,s in zip(Point,[x,y,z]):
            divExpr = divExpr.subs([(s,p)])
        return divExpr
def Curl(Field,Point=None,VarList=None):
    """
    Computes the curl of a vector field.

    The curl of a vector field F = <P, Q, R> in 3D is the vector field obtained by taking 
    the cross product of the del operator with the vector field.

    Parameters:
    Field (list of sympy expressions): The components of the vector field (e.g., [P, Q, R]).
    Point (list of sympy expressions, optional): A specific point at which to evaluate the curl.
    VarList (list of sympy symbols, optional): The variables with respect to which to compute the curl. Defaults to [x, y, z].

    Returns:
    list: The curl of the vector field, or its value at the specified point.

    Example:
    >>> Field = [x**2, y**2, z**2]
    >>> Curl(Field)
    [2*z - 2*y, 2*x - 2*z, 2*y - 2*x]
    """
    VarList = VarList or [x, y, z]
    gList = [[diff(f, v) for v in VarList] for f in Field]
    
    if len(gList) == 3:
        curl_array = [
            sp.trigsimp(gList[2][1] - gList[1][2]),
            sp.trigsimp(gList[0][2] - gList[2][0]),
            sp.trigsimp(gList[1][0] - gList[0][1])
        ]
    else:
        curl_array = [0, 0, sp.trigsimp(gList[1][0] - gList[0][1])]

    if Point is None:
        return curl_array
    else:
        curlEval = [(k.subs(i,j)) for i,j,k in zip(VarList,Point,curl_array)]
        return curlEval
def FindPotentialFnction(Field,VarList=None):
    if Curl(Field,VarList=VarList).tolist()!=[0]*len(Field):
        raise ValueError('Field:',Field,'Is not conservative. Therefore such a function does not exist.')
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    NewList = []
    for f,s in zip(Field,VarList):
        for i in sp.expand(sp.integrate(f,s)).as_ordered_terms():
            if i not in NewList:
                NewList.append(i)
    potenExpr = sum([i for i in NewList])
    return potenExpr
def ParametrizeExpr(Expression,MappedArray,VarList=None):
    """
    Substitutes variables in an expression using a mapped array or variable list.

    This function replaces the symbols in the expression with the corresponding values
    from the MappedArray or VarList.

    Parameters:
    Expression (sympy expression): The mathematical expression to be parametrized.
    MappedArray (list or sympy.Array): A list of values to substitute into the expression.
    VarList (list of sympy symbols, optional): The variables in the expression that are to be replaced. Defaults to None.

    Returns:
    sympy expression: The parametrized expression after substitution.

    Raises:
    ValueError: If MappedArray and VarList do not have the same length.

    Example:
    >>> Expression = x**2 + y**2
    >>> MappedArray = [1, 2]
    >>> ParametrizeExpr(Expression, MappedArray)
    5
    """
    if VarList is not None:
        if len(MappedArray) != len(VarList):
            raise ValueError('MappedArray and VarList must have the same length.')
        return Expression.subs(zip(VarList, MappedArray))
    return Expression.subs(zip(Expression.free_symbols, MappedArray))
def CountourInt(Field,VarList,VarBounds=None,Eval=False,Hard=True,ParaField=None,ParVar=None,ParVars=None,ParVarBounds=None,):
    """
    Computes the contour integral for a given vector field.

    This function calculates the contour integral of a vector field along a parametrized curve, 
    using Stokes' theorem or direct integration, depending on the input parameters.

    Parameters:
    Field (list of sympy expressions): The components of the vector field.
    VarList (list of sympy symbols): The variables in the vector field.
    VarBounds (list, optional): The bounds for the integration variables.
    Eval (bool, optional): Whether to evaluate the integral immediately. Defaults to False.
    Hard (bool, optional): Whether to use a parametrized approach or a simpler one. Defaults to True.
    ParaField (list, optional): The parametrized field for the curve.
    ParVar (sympy symbol, optional): The parameter variable used for parametrization.
    ParVars (list, optional): A list of parametrization variables for the curve.
    ParVarBounds (list, optional): The bounds for the parametrization variables.

    Returns:
    sympy expression: The contour integral expression or its evaluated result.

    Example:
    >>> Field = [x**2, y**2, z**2]
    >>> VarList = [x, y]
    >>> CountourInt(Field, VarList, Eval=True)
    Integral result
    """
    if Hard == True and ParVars == None:
        ParametrizedField = []
        for i in Field:
            ParametrizedField.append(i.subs(zip(VarList,ParaField)))
        integrand = 0
        for i in range(len(ParaField)):
            integrand = integrand + ParametrizedField[i]*sp.diff(ParaField[i],ParVar)
        hardIntegral = sp.Integral(integrand,(ParVar,ParVarBounds[0]))
        if Eval == False:
            print('Integrand for CounterInt =\n',integrand)
            return hardIntegral
        else:
            return hardIntegral.doit()
    else:
        if VarBounds==None:
            StokesIntegrand = (InteriorProduct(Curl(Field,VarList=VarList),[1,1,1])).subs(zip(VarList,ParaField))
            if ParaField == [sp.Symbol('r')*sp.cos(sp.Symbol('theta')),sp.Symbol('r')*sp.sin(sp.Symbol('theta'))] or ParaField == [sp.Symbol('r')*sp.cos(sp.Symbol('theta')),sp.Symbol('r')*sp.sin(sp.Symbol('theta')),sp.Symbol('z')]:
                StokesIntegrand *= sp.Symbol('r')
            StokesIntegral = sp.Integral(sp.simplify(StokesIntegrand),(ParVars[0],ParVarBounds[0]),(ParVars[1],ParVarBounds[1]))
            if Eval:
                return StokesIntegral.doit()
            print('Integrand for ContourInt =\n',StokesIntegrand)
            return StokesIntegral
        StokesIntegrand = InteriorProduct(Curl(Field,VarList=VarList),[1,1,1])
        StokesIntegral = sp.Integral(StokesIntegrand,(VarList[0],VarBounds[0]),(VarList[1],VarBounds[1]))
        if Eval == False:
            print('Integrand for ContourInt =\n',StokesIntegrand)
            return StokesIntegral
        return StokesIntegral.doit()
def SurfArea(Surface,VarList,ParaField=list|None,ParVarList=list|None,Eval=True,ParVarBounds=list|None,Easy=False,Left=True,IntExpr=None):
    """
    Computes the surface area of a parametrized surface.

    This function calculates the area of a surface using a double integral, either in an easy or full form.

    Parameters:
    Surface (list of sympy expressions): The components of the surface to be parametrized.
    VarList (list of sympy symbols): The variables used in the surface parametrization.
    ParaField (list, optional): The parametrized field for the surface.
    ParVarList (list, optional): The parametrization variables for the surface.
    Eval (bool, optional): Whether to evaluate the surface area immediately. Defaults to True.
    ParVarBounds (list, optional): The bounds for the parametrization variables.
    Easy (bool, optional): Whether to use an easy method for surface area computation. Defaults to False.
    Left (bool, optional): Determines if the left or right side of the surface is used for calculation.
    IntExpr (list, optional): The integration expression used for calculation.

    Returns:
    sympy expression: The surface area expression or its evaluated result.

    Example:
    >>> Surface = [x, y, z]
    >>> SurfArea(Surface, [x, y])
    Surface area result
    """
    if Easy == False:
        ParaSurf = []
        for i in Surface:
            ParaSurf.append(i.subs(zip(VarList,ParaField)))
        NormalPre = []
        for j in range(2):
            NormalPre.append([sp.trigsimp(diff(ParaSurf[i],ParVarList[j])) for i in range(3)])
        Integrand = sp.factor(Magnitude(CrossProduct(NormalPre[0],NormalPre[1])))
        VarList = ParVarList
    else:
        Integrand = sp.sqrt(diff(Surface[2],VarList[0])**2+diff(Surface[2],VarList[1])**2+1)
    if IntExpr != None:
        Integrand = Integrand*IntExpr
        if ParVarList != None:
            VarList=ParVarList
    if Eval == True:
        if Left == True:
            Integral = sp.Integral(Integrand,(VarList[0],ParVarBounds[0]),(VarList[1],ParVarBounds[1]))
        else:
            Integral = sp.Integral(Integrand,(VarList[1],ParVarBounds[1]),(VarList[0],ParVarBounds[0]))
        print('Surface Integral integrand is: \n',Integrand,'\nSurface Integral is:\n',Integral)
        EvalDIntegral = Integral.doit()
        return EvalDIntegral
    else:
        return Integrand
def SurfFlux(Surface,Field,VarList,ParaVarList,Bounds=list|None,Eval=True,Left=True,Upward=False):
    """
    Computes the surface flux of a vector field through a parametrized surface.

    The surface flux is the integral of the dot product between the vector field and the normal vector to the surface. 
    This function calculates the flux using the given parametrization of the surface and the vector field.

    Parameters:
    Surface (list of sympy expressions): The components of the surface, which are functions of the parametrization variables.
                                          For example, the surface might be given as [x(u,v), y(u,v), z(u,v)].
    Field (list of sympy expressions): The components of the vector field to compute the flux for. For example, [P(u,v), Q(u,v), R(u,v)].
    VarList (list of sympy symbols): The variables in the surface parametrization (e.g., [u, v]) used to describe the surface.
    ParaVarList (list of sympy symbols): The list of parametrization variables (e.g., [u, v]) that parameterize the surface.
    Bounds (list of tuples, optional): The bounds for the parametrization variables. If None, the function will compute flux over the entire domain.
    Eval (bool, optional): If True, the function will return the evaluated flux value. Defaults to True.
    Left (bool, optional): If True, the integration will be performed with respect to the first variable in ParaVarList as the outermost integral.
                           If False, the second variable is used as the outermost integral. Defaults to True.
    Upward (bool, optional): If True, the flux will be computed considering the upward direction of the normal vector (negative if downward).
                             Defaults to False, where the flux is computed for the default orientation.

    Returns:
    sympy expression or float: The computed surface flux, or the flux integrand if Eval=False.

    Raises:
    ValueError: If the dimensions of the surface and vector field do not match.

    Example:
    >>> Surface = [r*sp.cos(theta), r*sp.sin(theta), z]
    >>> Field = [x, y, z]
    >>> ParaVarList = [r, theta]
    >>> SurfFlux(Surface, Field, [x, y, z], ParaVarList, Eval=True)
    pi
    """
    ParaField = [f.subs(zip(VarList, Surface)) for f in Field]
    NormalDiffs = [[sp.trigsimp(diff(Surface[i], ParaVarList[j])) for i in range(3)] for j in range(2)] # gets the partial derivatives of the normal vector.
    Integrand = sp.factor(InteriorProduct(CrossProduct(NormalDiffs[0],NormalDiffs[1]),ParaField))
    if Upward:
        Integrand *= -1
    if Eval:
        if Left:
            Integral = sp.Integral(Integrand,(ParaVarList[0],Bounds[0]),(ParaVarList[1],Bounds[1]))
        else:
            Integral = sp.Integral(Integrand,(ParaVarList[1],Bounds[1]),(ParaVarList[0],Bounds[0]))
        print('Surface Integral integrand is: \n',Integrand,'\nSurface Integral is:\n',Integral)
        EvalDIntegral = Integral.doit()
        return EvalDIntegral
    else:
        return Integrand

#print(Jacobian(u*sp.cos(v),u*sp.sin(v)))
#matrix1 = sp.Matrix([[diff(u*sp.sin(v)*sp.cos(w),u),diff(u*sp.sin(v)*sp.cos(w),v)],[diff(u*sp.sin(v)*sp.sin(w),u),diff(u*sp.sin(v)*sp.sin(w),v)]])
#matrix1 = matrix1.col_insert(2, sp.Matrix([[diff(u*sp.sin(v)*sp.cos(w),w)],[diff(u*sp.sin(v)*sp.sin(w),w)]]))
#print(matrix1)
#Testing with the matrix created by the function.
#print(Jacobian(u*sp.sin(v)*sp.cos(w),u*sp.sin(v)*sp.sin(w),u*sp.cos(v)))
assert Jacobian([u*sp.cos(v),u*sp.sin(v)],[u,v]) == u , 'The jacobian of a polar transformation should equal r = u.'
assert Jacobian([u*sp.sin(v)*sp.cos(w),u*sp.sin(v)*sp.sin(w),u*sp.cos(v)],[u,v,w]) == (u**2)*sp.sin(v), 'The jacobian of the mapping from cartesian to spherical co-ordinates should equal p^2sin(phi) = u^2sin(v)'
