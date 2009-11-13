# -*- coding: cp1251 -*-
import math
def DataForm (AllDataArray):
    '''Назначение:
       формирование мачссива информативных параметров.
  Входные данные:
      SectMinCount - минимальное коичество сечений в испытаниях данного типа;
      NTest - количество испытаний;
  Выходные параметры:
      DataMCArray - массив информативных прищнаков;
      SectionNumbers - номера сечений, в которых берутся информативные признаки;
      HArray - энтропия;
      кол-во информативных параметров.
    '''
    SectMinCount = len(AllDataArray[0])
    NTest = len(AllDataArray)
    K = 5  #кол-во ближайших соседей для вычисления вероятности
    N = 0
    TMP = []
    for i in range(0, NTest):
        TMP.append(0.0)
    H1 = []
    for i in range(0, SectMinCount):
        H1.append(0.0)
    ParamsCount1 = len(AllDataArray[0][0])
    HArray=[]
    DataMCArray = []
    for i in range(0, NTest):
        DataMCArray.append([])
    for i in range(0, ParamsCount1):
        HArray.append([])
        for j in range(0, SectMinCount):
            sum = 0.0
            #формируем массив для i-го параметра и j-го сечения по всем испытаниям
            #for (int k = 0; k < NTest; k++)
            for k in range(0, NTest): 
                TMP[k] = AllDataArray[k][j][i]
            #сортируем массив по возрастанию
            for k in range(0, NTest-1):
                for m in range(k+1, NTest):
                    if (TMP[m] < TMP[k]):
                        c = TMP[k]
                        TMP[k] = TMP[m]
                        TMP[m] = c
            delta = NTest-1
            for k in range(0, NTest):
                Dist = []
                for cou in range(0, K):
                    Dist.append(0.0)
                m = 0
                s01 = 1
                s02 = 1
                sf = 0
                se = 0
                while (m < K): # может быть стоит здесь убрать <= и просто оставить <
                    tmp1 = 0.0
                    tmp2 = 0.0
                    s1 = s01
                    s2 = s02
                    #идем влево
                    while (((k - s1) >= 0) and not tmp1):
                        tmp = TMP[k] - TMP[k-s1]
                        if (tmp):
                            if(not m or ((m > 0) and (tmp > Dist[m-1]))):
                                tmp1 = tmp
                        s1 = s1 + 1

                    #идем вправо
                    while (((k + s2) < NTest) and  not tmp2):
                        tmp = TMP[k+s2] - TMP[k]
                        if (tmp):
                            if ((not m) or ((m > 0) and (tmp > Dist[m-1]))):
                                tmp2 = tmp
                        s2 = s2 + 1

                    if ((not tmp1) or (tmp1 and tmp2 and (tmp1 > tmp2))):
                        Dist[m] = tmp2
                        s02 = s2
                        se = s2-1

                    if ((not tmp2) or (tmp1 and tmp2 and (tmp2 > tmp1))):
                        Dist[m] = tmp1
                        s01 = s1
                        sf = s1-1

                    if (tmp1 and tmp2 and (tmp1 == tmp2)):
                        Dist[m] = tmp1
                        s01 = s1
                        s02 = s2
                        sf = s1-1
                        se = s2-1
                    m = m + 1

        
            rrr = 0
            #вычислим количество точек на отрезке [k-Dist; k+Dist]
            while ((TMP[k+se] - TMP[k]) >= (TMP[k] - TMP[k-sf])):
                sf = sf + 1
                rrr = 1
                if ((k - sf) < 0):
                    break

            if (rrr):
                sf = sf - 1
            rrr = 0
            while ((TMP[k+se] - TMP[k]) <= (TMP[k] - TMP[k-sf])):
                se = se + 1
                rrr = 1
                if ((k + se) >= NTest):
                    break

            if (rrr):
                se = se - 1
            delta_i = sf + se
            p = delta_i / float(delta)
            sum = sum + p * math.log10(p)
            HArray[i].append(-sum)
            H1[j] = -sum;


        nH = 3
        mH = [0, 0, 0]
        for j in range(0, nH):
            max = -1000
            for k in range(1, SectMinCount):
                if (H1[k] > max):
                    max = H1[k]
                    mH[j] = k
            H1[mH[j]] = -1000


        for k in range(0, nH):
            for j in range(0, NTest):
                DataMCArray[j].append(AllDataArray[j][mH[k]][i])
                #SectionNumbers[N] = mH[k];
            N = N + 1
    return DataMCArray
    
    
def MainComponents(DataArray, DimI, DimJ):
    '''Назначение:
       формирование мачссива координат в главных компонентах.
  Входные данные:
      DataArray - массив информативных прищнаков;
      DimI - число строк массива Х;
      DimJ - число столбцов массива Х;
  Выходные параметры:
      EV - собственные векторы;
      EN - собственные числа;
      Cord - координаты в главных компонентах.
    '''

  #
    DimI = len(DataArray)
    DimJ = len(DataArray[0])
    DataArray = Norm(DimI,DimJ,DataArray) #нормировка на среднеквадратичное отклонение
    CR = [] #корр. матрица
    for i in range(0, DimJ):
        CR.append([])
        for j in range(0, DimJ):
            CR[i].append(0.0)
    CR = Correl(DimI,DimJ,DataArray,CR); #вычисление корреляционной матрицы
    EV = Eigen_Vect(DimJ,CR) #определение собств. чисел и векторов
  #

    norm = []
    for i in range(0, DimJ):
        norm.append(0.0)
    for k in range(0, DimJ):
        for i in range(0, DimJ):
            norm[k] = norm[k] + EV[i][k] * EV[i][k]
    

    for i in range(0, DimJ):
        norm[i] = math.sqrt(norm[i])
    for k in range(0, DimJ):
        for i in range(0, DimJ):
            EV[i][k] = EV[i][k]/ float(norm[k])
  #проецирование данных в пространство собственных векторов
    Cord = []
    for i in range(0, DimI):
        Cord.append([])
        for j in range(0, DimJ):
            Cord[i].append(0.0)
        for k in range(0, DimJ):
            for j in range(0, DimJ):
                Cord[i][j] = Cord[i][j] + DataArray[i][k]*EV[k][j]
    return Cord

def Norm(DimI, DimJ, X):
    '''
  Назначение:
    нормирование выборки данных, содержащейся во входном массиве Х,
    на среднеквадратичное отклонение.
  Параметры:
    DimI - число строк массива Х;
    DimJ - число столбцов массива Х;
    '''
    #oпределение вектора средних
    for j in range(0, DimJ):
        MV = 0.
        D = 0.
        for i in range(0, DimI):
            MV = MV + X[i][j]
        MV = MV/float(DimI)
        #oпределение дисперсий
        for i in range(0, DimI):
            D = D + (X[i][j]-MV)*(X[i][j]-MV)
        D = math.sqrt(D)
        D = D / float(DimI)
        #cобственно нормирование
        for i in range(0, DimI):
            X[i][j] = (X[i][j]-MV)/float(D)

    return X

def Correl (DimI, DimJ, X, CR):
    '''
  Назначение:
    вычисление корреляционной матрицы для
    выборки данных, содержащейся во входном массиве Х.
  Параметры:
    DimI - число строк массива Х;
    DimJ - число столбцов массива Х;
    CR - корреляционная матрица [DimJ][DimJ];
    MV - вектор средних;
    D - вектор дисперсий.
    '''
 
    MV = []
    D = []
    for i in range(0, DimJ):
        MV.append(0.0)
        D.append(0.0)

    TMP = []
    for i in range(0, DimI):
        TMP.append([])
        for j in range(0, DimJ):
            TMP[i].append(0.0)
            
    for i in range(0, DimI):
        for j in range(0, DimJ):
            TMP[i][j] = X[i][j]

    for j in range(0, DimJ):
        MV[j] = 0.0

        for i in range(0, DimI):
            MV[j] = MV[j] + TMP[i][j]
        MV[j] = MV[j]/float(DimI)
        D[j] = 0.0

        for i in range(0, DimI):
            D[j] = D[j]+ (TMP[i][j] - MV[j]) * (TMP[i][j] - MV[j])
        D[j] = D[j]/float(DimI)
        D[j] = math.sqrt(D[j])
        for i in range(0, DimI):
            TMP[i][j] = TMP[i][j] - MV[j]

    for i in range(0, DimJ-1):
        for j in range(i+1, DimJ):
            CR[i][j] = 0.0
            for k in range(0, DimI):
                CR[i][j] = CR[i][j] + TMP[k][j] * TMP[k][i]
            CR[i][j] = CR[i][j]/float(DimI)
            CR[i][j] = CR[i][j]/float(D[i] * D[j])
            CR[j][i] = CR[i][j]


    for i in range(0, DimJ):
        CR[i][i] = 1.0

    return CR
    
    
def Eigen_Vect(DimJ, CR):

    '''
  Назначение:
    вычисление собственных значений и собственных векторов
    симметричной квадратной матрицы CR методом вращения.
  Параметры:
    Dim - размерность массива CR;
    EV - матрица [DimJ][DimJ] компонент собственных векторов;
    EN - собственные значения.
    '''
    R = 0.1E-6
    S = []
    F = []
    EN = []
    EV = []
    for i in range(0, DimJ):
        S.append([])
        EV.append([])
        F.append(0)
        EN.append(0.0)
        for j in range(0, DimJ):
            S[i].append(0.0)
            EV[i].append(0.0)

    for i in range(0, DimJ):
        for j in range(0, DimJ):
            S[i][j] = 0

    for i in range(0, DimJ):
        S[i][i] = 1

    c1 = 0.0
    c2 = 0.0
    n1 = 0.0
    n2 = 0.0
    t = 0.0
    t1 = 0.0
    t2 = 0.0
    t3 = 0.0
    v1 = 0.0
    v2 = 0.0
    v3 = 0.0
    m1 = 0.0
    i1 = 0.0
    w = 0.0
    for i in range(1, DimJ):
        for j in range(0, i):
            i1 += 2*CR[i][j]*CR[i][j]
    n1 = math.sqrt(i1)
    n2 = (R/float(DimJ))*n1
    t = n1
    i2 = 1
    while (t > n2):
        t = t/float(DimJ)
        while (i2):
            i2 = 0
            for q in range(1, DimJ):
                for p in range(0, q):
                    if (math.fabs(CR[p][q]) <= t):
                        continue
                    i2 = 1
                    v1 = CR[p][p]
                    v2 = CR[p][q]
                    v3 = CR[q][q]
                    m1 = (v1-v3)*0.5
                    if (m1):
                        w = (-sign(m1))*v2/math.sqrt(v2*v2+m1*m1)
                    else:
                        w = -1
                    t1 = w/math.sqrt(2*(1+math.sqrt(1-w/2)))
                    t2 = t1*t1
                    print t2
                    c1 = math.sqrt(1-t2)
                    c2 = c1*c1
                    t3 = t1*c1
                    for i in range(0, DimJ):
                        i1 = CR[i][p]*c1 - CR[i][q]*t1
                        CR[i][q] = CR[i][p]*t1 + CR[i][q]*c1
                        CR[i][p] = i1
                        i1 = S[i][p]*c1 - S[i][q]*t1
                        S[i][q] = S[i][p]*t1 + S[i][q]*c1
                        S[i][p] = i1

                    for i in range(0, DimJ):
                        CR[p][i] = CR[i][p]
                        CR[q][i] = CR[i][q]

                    CR[p][p] = v1*c2 + v3*t2 - 2*v2*t3
                    CR[q][q] = v1*t2 + v3*c2 + 2*v2*t3
                    CR[p][q] = (v1-v3)*t3 + v2*(c2-t2)
                    CR[q][p] = CR[p][q]


    for i in range(0, DimJ):
        EN[i] = CR[i][i]
        F[i] = i

    p = DimJ
    while (p):
        p = p - 1
        for i in range(0, p):
            if (EN[i] > EN[i+1]):
                continue
            c1 = EN[i]
            EN[i] = EN[i+1]
            EN[i+1] = c1
            j = F[i]
            F[i] = F[i+1]
            F[i+1] = j


    for i in range(0, DimJ):
        k = F[i]
        for j in range(0, DimJ):
            EV[j][i] = S[j][k]
    return EV

def sign(X):
    '''
  Назначение:
    вычисление знака числа Х.
    '''

    if (math.fabs(X) < 1.0e-9):
        return 0
    else:
        if (X < 0):
            return -1
        else:
            return 1


