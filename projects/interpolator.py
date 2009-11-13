# -*- coding: utf-8 -*-
def buildlinearspline(x, y, n):
    #
    # Sort points
    #
    heapsortpoints(x, y, n);
    
    #
    # Fill C:
    #  C[0]            -   length(C)
    #  C[1]            -   type(C):
    #                      3 - general cubic spline
    #  C[2]            -   N
    #  C[3]...C[3+N-1] -   x[i], i = 0...N-1
    #  C[3+N]...C[3+N+(N-1)*4-1] - coefficients table
    #
    tblsize = int(3+n+(n-1)*4)
    c = []
    for i in range(0, tblsize):
        c.append(0.0)
    c[0] = tblsize
    c[1] = 3
    c[2] = n
    for i in range(0, n):
        c[3+i] = x[i]
    for i in range(0, n-1):
        c[3+n+4*i+0] = y[i]
        try:
            c[3+n+4*i+1] = (y[i+1]-y[i])/float((x[i+1]-x[i]))
        except:
            c[3+n+4*i+1] = (y[i+1]-y[i])/float((x[i+1]-x[i]+1))
        c[3+n+4*i+2] = 0
        c[3+n+4*i+3] = 0
    return c
'''
Построение таблицы коэффициентов сплайна Эрмита

Входные параметры:
    X           -   абсциссы, массив с нумерацией элементов [0..N-1]
    Y           -   значения функции,
                    массив с нумерацией элементов [0..N-1]
    D           -   значения производной,
                    массив с нумерацией элементов [0..N-1]
    N           -   число точек, N>=2

Выходные параметры:
    C           -   таблица коэффициентов сплайна для использования в
                    подпрограмме SplineInterpolation

'''
def buildhermitespline(x, y, d, n):

    if (n < 2):
        raise Exception('The number of points is not enought!')
    
    #
    # Sort points
    #
    heapsortdpoints(x, y, d, n);
    
    #
    # Fill C:
    #  C[0]            -   length(C)
    #  C[1]            -   type(C):
    #                      3 - general cubic spline
    #  C[2]            -   N
    #  C[3]...C[3+N-1] -   x[i], i = 0...N-1
    #  C[3+N]...C[3+N+(N-1)*4-1] - coefficients table
    #
    c = []
    tblsize = 3+n+(n-1)*4;
    for i in range(0, tblsize):
        c.append(0)
    c[0] = tblsize;
    c[1] = 3;
    c[2] = n;
    for i in range(0, n):
        c[3+i] = x[i];

    for i in range(0, n-1):
        delta = x[i+1]-x[i];
        delta2 = math.sqrt(delta);
        delta3 = delta*delta2;
        c[3+n+4*i+0] = y[i];
        c[3+n+4*i+1] = d[i];
        c[3+n+4*i+2] = (3*(y[i+1]-y[i])-2*d[i]*delta-d[i+1]*delta)/float(delta2);
        c[3+n+4*i+3] = (2*(y[i]-y[i+1])+d[i]*delta+d[i+1]*delta)/float(delta3);

    return c


'''
Построение таблицы коэффициентов сплайна Акимы

Входные параметры:
    X           -   абсциссы, массив с нумерацией элементов [0..N-1]
    Y           -   значения функции,
                    массив с нумерацией элементов [0..N-1]
    N           -   число точек, N>=5

Выходные параметры:
    C           -   таблица коэффициентов сплайна для использования в
                    подпрограмме SplineInterpolation

'''
import math
def buildakimaspline(x, y, n):
    
    d = []
    w = []
    diff = []
    for i in range(0, n):
        d.append(0)
        w.append(0)
        diff.append(0)

    if (n < 5):
        raise Exception('The number of points is not enought!')
    

    # Sort points
    heapsortpoints(x, y, n)
    
    #
    # Prepare W (weights), Diff (divided differences)
    #
    for i in range(0, n-1):
        diff[i] = (y[i+1]-y[i])/float((x[i+1]-x[i]))

    for i in range(1, n-1):
        w[i] = math.fabs(diff[i]-diff[i-1])
    
    #
    # Prepare Hermite interpolation scheme
    #
    for i in range(2, n-2):
        if( (math.fabs(w[i-1])+math.fabs(w[i+1]))!=0 ):
            d[i] = (w[i+1]*diff[i-1]+w[i-1]*diff[i])/float((w[i+1]+w[i-1]))
        else:
            d[i] = ((x[i+1]-x[i])*diff[i-1]+(x[i]-x[i-1])*diff[i])/float((x[i+1]-x[i-1]))


    d[0] = diffthreepoint(x[0], x[0], y[0], x[1], y[1], x[2], y[2])
    d[1] = diffthreepoint(x[1], x[0], y[0], x[1], y[1], x[2], y[2])
    d[n-2] = diffthreepoint(x[n-2], x[n-3], y[n-3], x[n-2], y[n-2], x[n-1], y[n-1])
    d[n-1] = diffthreepoint(x[n-1], x[n-3], y[n-3], x[n-2], y[n-2], x[n-1], y[n-1])

    
    #
    # Build Akima spline using Hermite interpolation scheme
    #
    res = buildhermitespline(x, y, d, n)
    return res


'''
Вычисление интерполирующего сплайна

Входные параметры:
    C   -   массив коэффициентов, вычисленный подпрограммой для
            построения сплайна.
    X   -   точка, в которой вычисляется значение сплайна

Результат:
    значение сплайна в точке X

'''

def splineinterpolation(c, x):

    if (round(c[1])!=3):
        raise Exception ("SplineInterpolation: incorrect C!")
    n = round(c[2])
    
    #
    # Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    #
    l = 3
    r = 3+n-2+1
    while(l!=r-1):
        m = int((l+r)/2)
        if( c[m]>=x ):
            r = m
        else:
            l = m
    
    #
    # Interpolation
    #
    x = x-c[l]
    m = int(3+n+4*(l-3))
    result = c[m]+x*(c[m+1]+x*(c[m+2]+x*c[m+3]))
    return result


'''
Internal subroutine. Heap sort.
'''
def heapsortpoints(x, y, n):
    
    #
    # Test for already sorted set
    #
    isascending = True;
    isdescending = True;

    for i in range(1, n-1):
        isascending = isascending and (x[i] > x[i-1])
        isdescending = isdescending and (x[i] < x[i-1])

    if( isascending ):
        return True

    if( isdescending ):
        for i in range(0, n-1):
            j = n-1-i
            if( j<=i ):
                break
            tmp = x[i]
            x[i] = x[j]
            x[j] = tmp
            tmp = y[i]
            y[i] = y[j]
            y[j] = tmp
        return True
    
    #
    # Special case: N=1
    #
    if(n == 1):
        return True


    #
    # General case
    #
    i = 2
    while(i<=n):
        t = i
        while(t != 1):
            k = t/2
            if( x[k-1] >= x[t-1] ):
                t = 1
            else:
                tmp = x[k-1]
                x[k-1] = x[t-1]
                x[t-1] = tmp
                tmp = y[k-1]
                y[k-1] = y[t-1]
                y[t-1] = tmp
                t = k
        i = i+1

    i = n-1
    while(i>=1):
        tmp = x[i]
        x[i] = x[0]
        x[0] = tmp
        tmp = y[i]
        y[i] = y[0]
        y[0] = tmp
        t = 1
        while(t!=0):
            k = 2*t
            if( k>i ):
                t = 0
            else:
                if( k<i ):
                    if( x[k]>x[k-1] ):
                        k = k+1
                if( x[t-1]>=x[k-1] ):
                    t = 0
                else:
                    tmp = x[k-1]
                    x[k-1] = x[t-1]
                    x[t-1] = tmp
                    tmp = y[k-1]
                    y[k-1] = y[t-1]
                    y[t-1] = tmp
                    t = k
        i = i-1

        
'''
Internal subroutine. Heap sort.
'''
def heapsortdpoints(x, y, d, n):
    
    #
    # Test for already sorted set
    #
    isascending = True
    isdescending = True
    for i in range(1, n):
        isascending = isascending and (x[i]>x[i-1])
        isdescending = isdescending and (x[i]<x[i-1])

    if( isascending ):
        return True

    if( isdescending ):
        for i in range(0, n):
            j = n-1-i
            if( j<=i ):
                break
            tmp = x[i]
            x[i] = x[j]
            x[j] = tmp
            tmp = y[i]
            y[i] = y[j]
            y[j] = tmp
            tmp = d[i]
            d[i] = d[j]
            d[j] = tmp
        return True
    
    #
    # Special case: N=1
    #
    if( n==1 ):
        return True
    
    #
    # General case
    #
    i = 2
    while(i<=n):
        t = i
        while(t!=1):
            k = t/2
            if( x[k-1]>=x[t-1] ):
                t = 1
            else:
                tmp = x[k-1]
                x[k-1] = x[t-1]
                x[t-1] = tmp
                tmp = y[k-1]
                y[k-1] = y[t-1]
                y[t-1] = tmp
                tmp = d[k-1]
                d[k-1] = d[t-1]
                d[t-1] = tmp
                t = k
        i = i+1

    i = n-1;
    while(i>=1):
        tmp = x[i]
        x[i] = x[0]
        x[0] = tmp
        tmp = y[i]
        y[i] = y[0]
        y[0] = tmp
        tmp = d[i]
        d[i] = d[0]
        d[0] = tmp
        t = 1
        while(t!=0):
            k = 2*t
            if( k>i ):
                t = 0
            else:
                if( k<i ):
                    if( x(k)>x(k-1) ):
                        k = k+1
                if( x(t-1)>=x(k-1) ):
                    t = 0
                else:
                    tmp = x[k-1]
                    x[k-1] = x[t-1]
                    x[t-1] = tmp
                    tmp = y[k-1]
                    y[k-1] = y[t-1]
                    y[t-1] = tmp
                    tmp = d[k-1]
                    d[k-1] = d[t-1]
                    d[t-1] = tmp
                    t = k
        i = i-1


'''
Internal subroutine. Three-point differentiation
'''
def diffthreepoint(t, x0, f0, x1, f1, x2, f2):

    t = t-x0
    x1 = x1-x0
    x2 = x2-x0
    a = (f2-f0-x2/float(x1*(f1-f0)))/float((math.sqrt(x2)-x1*x2))
    b = (f1-f0-a*math.sqrt(x1))/float(x1)
    result = 2*a*t+b;
    return result
    
'''
Распаковка сплайна

Входные параметры:
    C   -   массив коэффициентов, вычисленный подпрограммой для
            построения сплайна.

Выходные параметры:
    N   -   число точек, на основе которых был построен сплайн
    Tbl -   таблица коэффициентов сплайна. Массив с нумерацией элементов
            [0..N-2, 0..5].
            Для I = 0..N-2:
                Tbl[I,0] = X[i]
                Tbl[I,1] = X[i+1]
                Tbl[I,2] = C0
                Tbl[I,3] = C1
                Tbl[I,4] = C2
                Tbl[I,5] = C3
            Сплайн имеет вид:
                t = x-x[i]
                S(x) = C0 + C1*t + C2*t^2 + C3*t^3

'''
def splineunpack(c):
    n = round(c[2])
    tbl = []
    for i in range(0, n-1):
        tbl.append(0)
        for j in range(0, 6):
            tbl[i] = [0, 0, 0, 0, 0, 0] 
            
    
    #
    # Fill
    #
    for i in range(0, n-1):
        tbl[i][0] = c[3+i]
        tbl[i][1] = c[3+i+1]
        tbl[i][2] = c[int(3+n+4*i)]
        tbl[i][3] = c[int(3+n+4*i+1)]
        tbl[i][4] = c[int(3+n+4*i+2)]
        tbl[i][5] = c[int(3+n+4*i+3)]
        
    return tbl
    
 

def solveTriDiag(TDM, F, n, b):
    alph = []
    beta = []
    for i in range(0, n-1):
        alph.append(0.0)
        beta.append(0.0)
 
    
    alph[0] = - (TDM[2][0]/float(TDM[1][0]))
    beta[0] = (F[0]/float(TDM[1][0]))
 
    for i in range(1, n-1):
        alph[i] = -(TDM[2][i]/float((TDM[1][i] + TDM[0][i]*alph[i-1])))
        beta[i] = (F[i]-TDM[0][i]*beta[i-1])/float((TDM[1][i] + TDM[0][i]*alph[i-1]))

    b[n-1] = (F[n-1]-TDM[0][n-1]*beta[n-2])/float((TDM[1][n-1] + TDM[0][n-1]*alph[n-2]))
 
    for i in range(n-2, -1, -1):
         b[i] = b[i+1] * alph[i] + beta[i]

 
def buildSpline(KnotArray, n):

    a = []
    c = []
    d = []
    delta = []
    h = []
    TriDiagMatrix = [[], [], []]
    b = []
    f = []
    for i in range(0, n-1):
        a.append(0.0)
        c.append(0.0)
        d.append(0.0)
        delta.append(0.0)
        h.append(0.0)
    for i in range(0, n):
        b.append(0.0)
        f.append(0.0)
    for i in range(0, 3):
        for j in range(0, n):
            TriDiagMatrix[i].append(0.0)


 
    if (n<3):
        return -1
 
    x3 = KnotArray[2][0] - KnotArray[0][0]
    xn = KnotArray[n-1][0] - KnotArray[n-3][0]
 
    for i in range(0, n-1):
        a[i] = KnotArray[i][1]
        h[i] = KnotArray[i+1][0] - KnotArray[i][0]
        delta[i] = (KnotArray[i+1][1] - KnotArray[i][1])/float(h[i])
        if i > 0:
            TriDiagMatrix[0][i] = h[i]
            f[i] = 3*(h[i]*delta[i-1] + h[i-1]*delta[i])
        else:
            TriDiagMatrix[0][i] = x3
            f[i] = 0
    TriDiagMatrix[1][0] = h[0]
    TriDiagMatrix[2][0] = h[0]
    for i in range(1, n-1):
        TriDiagMatrix[1][i] = 2*(h[i] + h[i-1])
        TriDiagMatrix[2][i] = h[i]

    TriDiagMatrix[1][n-1] = h[n-2]
    TriDiagMatrix[2][n-1] = xn
    TriDiagMatrix[0][n-1] = h[n-2]

    i = n-1
    f[0] = ((h[0]+2*x3)*h[1]*delta[0] + (h[0]**2)*delta[1])/float(x3)
    f[n-1]=((h[i-1]**2)*delta[i-2]+(2*xn+h[i-1])*h[i-2]*delta[i-1])/float(xn)
     
    solveTriDiag(TriDiagMatrix,f, n, b)
 
    Coef = []
    for i in range(0, n-1):
        Coef.append([0,0,0,0])
    for j in range(0, n-1):
        d[j] = (b[j+1]+b[j]-2*delta[j])/float((h[j]*h[j]))
        c[j] = 2*(delta[j]-b[j])/h[j]-(b[j+1]-delta[j])/float(h[j])
        Coef[j][0] = a[j]
        Coef[j][1] = b[j]
        Coef[j][2] = c[j]
        Coef[j][3] = d[j]
    return Coef

def interpolate(KnotArray, Coef, x):
    i = 0
    while (KnotArray[i][0] < x):
        i = i+1
    i = i - 1
    return Coef[i][0] + Coef[i][1]*(x-KnotArray[i][0]) + Coef[i][2]*((x-KnotArray[i][0])**2) + Coef[i][3]*((x-KnotArray[i][0])**3)  
 

 


