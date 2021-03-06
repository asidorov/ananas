import interpolator
import math
import re
import processor
def test():
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [1, 5, 3, 7, 6, 5, 4, 3, 2, 1]
    n = 10
    c = interpolator.buildakimaspline(x, y, n)
    res = []
    for i in [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]:
        res.append(interpolator.splineinterpolation(c, i))
    
    res = interpolator.splineunpack(c)
    print res
    
def test1():
    xy = [[1,1], [2,5], [3,3], [4,7], [5,6], [6,5], [7,4], [8,3], [9,2], [10,1]]
    n = len(xy)
    c = interpolator.buildSpline(xy, n)
    res = []
    for i in [1.2, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]:
        res.append(interpolator.interpolate(xy, c, i))
    print res

def test2():
    f = open('D:\\LEARN\\AnaNas-Python\\1.txt')
    xy = []
    i = 0
    for line in f:
        data = re.findall('\d+\.*\d*\.*\d*', str(line))
        time_com = data[0].split('.')
        time = float(time_com[0])*3600 + float(time_com[1])*60 + float(time_com[2])
        xy.append([time, float(data[2])])
    int_time = xy[0][0]
    for i in xy:
        i[0] = i[0]-int_time
    n = len(xy)
    c = interpolator.buildSpline(xy, n)
    res = []
    for i in range(0, int(xy[int(len(xy)-1)][0])):
        res.append(interpolator.interpolate(xy, c, i))
    res_file = open('D:\\LEARN\\AnaNas-Python\\res.txt', 'w')
    for i in range(0, int(xy[int(len(xy)-1)][0])):
        res_file.write(str(i)+' '+str(res[i])+'\r\n')
    init_file = open('D:\\LEARN\\AnaNas-Python\\init.txt', 'w')
    for i in range(0,int(len(xy)-1)):
        init_file.write(str(xy[i][0])+' '+str(xy[i][1])+'\r\n')
        
def test3():
    f = open('D:\\LEARN\\AnaNas-Python\\1.txt')
    x = []
    y = []
    i = 0
    for line in f:
        data = re.findall('\d+\.*\d*\.*\d*', str(line))
        time_com = data[0].split('.')
        time = float(time_com[0])*3600 + float(time_com[1])*60 + float(time_com[2])
        x.append(time)
        y.append(float(data[2]))
    int_time = x[0]
    for i in range(0,len(x)):
        x[i] = x[i]-int_time
    
    n = len(x)
    c = interpolator.buildlinearspline(x, y, n)
    res = []
    for i in range(0, int(x[int(len(x)-1)])):
        res.append(interpolator.splineinterpolation(c, i))
    res_file = open('D:\\LEARN\\AnaNas-Python\\res.txt', 'w')
    print x[int(len(x)-1)]
    for i in range(0, int(x[int(len(x)-1)])):
        res_file.write(str(i)+' '+str(res[i])+'\r\n')
    init_file = open('D:\\LEARN\\AnaNas-Python\\init.txt', 'w')
    for i in range(0,int(len(x))):
        init_file.write(str(x[i])+' '+str(y[i])+'\r\n')
        
def test4():
    fileBase = '/home/alexey/aspirant/ananas/projects/1/'
    AllDataArray = []
    testNumbers = 0
    for file in range(1, 9):
        f = open(fileBase+str(file)+'.TXT')
        testMeasures = 0
        AllDataArray.append([])
        for line in f:
            data = re.findall('\d+\.*\d*\.*\d*', str(line))
            time_com = data[0].split('.')
            time = float(time_com[0])*3600 + float(time_com[1])*60 + float(time_com[2])
            if testMeasures == 0:
                int_time = time
            time = time - int_time
            data[0] = time
            AllDataArray[testNumbers].append(data)
            testMeasures = testMeasures + 1
        testNumbers = testNumbers + 1
        testMeasures = 0
    min = AllDataArray[0][len(AllDataArray[0])-1][0]
    max = min
    for i in range(0, len(AllDataArray)):
        if  (AllDataArray[i][len(AllDataArray[i])-1][0] < min):
            min = AllDataArray[i][len(AllDataArray[i])-1][0]
        if  (AllDataArray[i][len(AllDataArray[i])-1][0] > max):
            max = AllDataArray[i][len(AllDataArray[i])-1][0]
    res_file = open('/home/alexey/aspirant/ananas/res.TXT', 'w')
    res_file.write(str(AllDataArray[3]))
    res_file.close()
    numTimeSection = 30
    InterpolatedDataArray = []
    time_step = min/numTimeSection
    time = 0
    # initialize the InterpolatedDataArray
    for i in range(0, len(AllDataArray)):
        InterpolatedDataArray.append([])
        for j in range(0, numTimeSection):
            InterpolatedDataArray[i].append([])
            for k in range(0, len(AllDataArray[i][j])-1):
                InterpolatedDataArray[i][j].append(0.0)

    for i in range(0, len(AllDataArray)):
        for j in range(1, len(AllDataArray[i][0])):
            x = []
            y = []
            for k in range(0, len(AllDataArray[i])):
                x.append(float(AllDataArray[i][k][0]))
                y.append(float(AllDataArray[i][k][j]))
            n = len(x)
            print('i is '+str(i))
            try:
                c = interpolator.buildlinearspline(x, y, n)
            except:
                print x, y
                raise Exception('blaa')
            for k in range(0, numTimeSection):
                InterpolatedDataArray[i][k][j-1] = interpolator.splineinterpolation(c, time)
                time = time + time_step
            time = 0
    DataMCArray = processor.DataForm(InterpolatedDataArray)
    EigenNumbers = []
    EigenVectors = []
    for i in range(0, len(DataMCArray)):
        EigenNumbers.append(0.0)
        EigenVectors.append([])
        for j in range(0, len(DataMCArray)):
            EigenVectors[i].append(0.0)

    Cord = processor.MainComponents(DataMCArray,len(DataMCArray),len(DataMCArray[0]))
    res_file = open('/home/alexey/aspirant/ananas/res.TXT', 'w')
    for i in Cord:    
        res_file.write(str(i[0])+' '+str(i[1])+'\n')
    return Cord
test4()
