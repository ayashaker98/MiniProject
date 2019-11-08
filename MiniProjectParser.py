
from liberty.types import Group, select_cell, select_pin, select_timing_table
from liberty.parser import parse_liberty
from collections import namedtuple
import array
import numpy
from scipy.spatial.distance import euclidean
from scipy.interpolate import interp1d, interp2d
liberty_file = "/Users/noha/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())
fname = "/Users/noha/Desktop/rca4.rtlnopwr.v"

import re

maxFanout = 3 #input("Please Enter the maximum Fanout for any cell : ") 
maxDelay = 0.01 #input("Please Enter the maximum Delay for any cell : ") 

def findMiddle(input_list):
    return(int((len(input_list) - 1)/2))

moduledefintion = ""
inn = []
out = []
wire = []
strmessage = ""
functions = []
for code in open(fname,"r"):
    line = code
    if("input" in line):
        line = re.sub("input", '', line) 
        inn.append(line)
    elif("output" in line):
        line = re.sub("output", '', line) 
        out.append(line)
    elif("wire" in line):
        line = re.sub("wire", '', line) 
        wire.append(line)
    elif("endmodule" in line):
        strmessage = "end of module"
    elif ("module" in line):
        moduledefintion = line
    elif(line == "\n"):
        donothing = ""
    else:
        functions.append(line)


NumberofCells = len(functions)


# In[483]:


dictList = []
for code in functions:
    n = code.split(maxsplit=2)
    del n[0]
    dictList.append(n)
    


# In[484]:


for i in range(len(dictList)):
    for j in range(2):
        n =  re.sub("\n", '', dictList[i][j])
        dictList[i][j] = n
        n =  re.sub(";", '', dictList[i][j])
        dictList[i][j] = n

for i in range(len(dictList)):
    n = dictList[i][1].split()
    dictList[i][1] = n
for i in range(len(dictList)):
    del dictList[i][1][0]
    del dictList[i][1][len(dictList[i][1])-1]

# print(dictList[0])
# print(dictList[0][0])
# print(NumberofCells)


start = '.'
end = '('
end2 = ')'
pins = namedtuple('pins', ['cell','pin', 'parameter', 'direction'])
allpins = []

for x in range(len(dictList)):
    for y in range(len(dictList[x])):
        for z in range(len(dictList[x][y])):
            if (dictList[x][y][z][dictList[x][y][z].find(start)+len(start):dictList[x][y][z].rfind(end)] != ""):
                allpins.append(pins(dictList[x][0], dictList[x][y][z][dictList[x][y][z].find(start)+len(start):dictList[x][y][z].rfind(end)], dictList[x][y][z][dictList[x][y][z].find(end)+len(end):dictList[x][y][z].rfind(end2)], "none"))


for f in range(len(allpins)-1):
    if (allpins[f].cell == allpins[f+1].cell):
        allpins[f] = allpins[f]._replace(direction = "input")
    else:
        allpins[f] = allpins[f]._replace(direction = "output")

allpins[len(allpins)-1] = allpins[len(allpins)-1]._replace(direction = "output")
# print(allpins)

Fanoutnumber = -1
count  = 0
FanoutCouple = namedtuple('Fanout', ['cell', 'fanout', 'status'])
Fanout = []
for l in range(len(allpins)):
    if (allpins[l].direction=="output"):
        count = count + 1
        for h in range(len(allpins)):
            if (allpins[h].parameter == allpins[l].parameter):
                Fanoutnumber = Fanoutnumber + 1
        Fanout.append(FanoutCouple(allpins[l].cell,Fanoutnumber,"none"))
        Fanoutnumber = -1

for g in range(len(Fanout)):
    if (Fanout[g].fanout > int(maxFanout)):
        Fanout[g] = Fanout[g]._replace(status = "Violating")
    else:
        Fanout[g] = Fanout[g]._replace(status = "Not Violating")


finalpins = namedtuple('finalpins', ['cellname', 'capacitance','pinwithmaxcap'])
CapPins = []
MaxCapPins = []
for i in range(len(Fanout)):
    for j in range(len(allpins)):
        if(Fanout[i].cell == allpins[j].cell):
            cgroup = library.get_group('cell', allpins[j].cell.split('_')[0])
            pingroup = cgroup.get_group('pin', allpins[j].pin)
            cap = pingroup.__getitem__('capacitance')
            CapPins.append(finalpins(allpins[j].cell, cap, allpins[j].pin))
    MaxCapPins.append(finalpins(CapPins[0].cellname, max(CapPins, key=lambda k: k.capacitance).capacitance, max(CapPins, key=lambda k: k.capacitance).pinwithmaxcap))
   # CapPins.clear()

import collections
Fanout2 = collections.namedtuple('Fanoutpairs', 'gate1 gate2 outputpin1 inputpin2 pin1 pin2 inputcap')
Fanoutpairs = []
for i in range(len(allpins)):
    for j in range(len(allpins)):
        if ((allpins[i].parameter  == allpins[j].parameter)&(allpins[i].direction == 'output')&(allpins[j].direction =='input')):
            # print("Parameter: " , allpins[i].parameter)
            # print("Gate 1 : ",allpins[i].cell)
            # print("Gate 2 : ",allpins[j].cell)
            Fanoutpairs.append(Fanout2(allpins[i].cell,allpins[j].cell,allpins[i].parameter,allpins[j].parameter,allpins[i].pin,allpins[j].pin,0))

for i in range(len(Fanoutpairs)):
    for j in range(len(CapPins)):
        if((CapPins[j].cellname == Fanoutpairs[i].gate2) & (CapPins[j].pinwithmaxcap == Fanoutpairs[i].pin2)):
            Fanoutpairs[i]= Fanoutpairs[i]._replace(inputcap = CapPins[j].capacitance)

FanoutSum = 0
Fanout4 = namedtuple('Fanout4', ['cell', 'load', 'outputpin', 'relatedpin'])
CellwithLoad = []
for f in range(len(Fanout)):
    for r in range(len(Fanoutpairs)):
        if (Fanout[f].cell == Fanoutpairs[r].gate1):
            FanoutSum = FanoutSum + Fanoutpairs[r].inputcap
    CellwithLoad.append(Fanout4(Fanout[f].cell, FanoutSum, 'none', 'none'))
    FanoutSum = 0
# print(CellwithLoad[0].load)




# print(findMiddle(timing_table_transition[0]))
# print(timing_table_capacitance)
# print(timing_table_values)

for u in range(len(CellwithLoad)):
    for g in range(len(allpins)):
        if(CellwithLoad[u].cell == allpins[g].cell):
            if(allpins[g].direction == "output"):
                CellwithLoad[u] = CellwithLoad[u]._replace(outputpin = allpins[g].pin)
            elif((allpins[g].direction == "input") and (CellwithLoad[u].relatedpin == 'none')):
                CellwithLoad[u] = CellwithLoad[u]._replace(relatedpin = allpins[g].pin)


rise_delay = -1
fall_delay = -1
minIndex1 = 0
minIndex2 = 0
distance = []
cell_delay_things = namedtuple('cell_delay_things', ['cell', 'delay'])
cell_delay = []

for i in range(len(CellwithLoad)):
    cellgroup = library.get_group('cell', CellwithLoad[i].cell.split('_')[0])
    pin = cellgroup.get_group('pin', CellwithLoad[i].outputpin)
    #for the rise delays
    timing_table_rise = select_timing_table(pin, CellwithLoad[i].relatedpin, 'cell_rise')
    timing_table_transition_rise = timing_table_rise.get_array("index_1")
    timing_table_capacitance_rise = timing_table_rise.get_array("index_2")
    timing_table_values_rise = timing_table_rise.get_array("values")
    #for the fall delays
    timing_table_fall = select_timing_table(pin, CellwithLoad[i].relatedpin, 'cell_fall')
    timing_table_transition_fall = timing_table_fall.get_array("index_1")
    timing_table_capacitance_fall = timing_table_fall.get_array("index_2")
    timing_table_values_fall = timing_table_fall.get_array("values")
    for j in range(len(timing_table_capacitance_rise[0])):
        if (CellwithLoad[i].load == timing_table_capacitance_rise[0][j]):
            rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][j]
            print("Rise Delay was set successfuly")
    if(rise_delay == -1):
        # print("Rise Delay not set, must do interpolation or extrapolation!!!")
        for k in range(len(timing_table_capacitance_rise[0])):
            distance.append(euclidean(CellwithLoad[i].load,timing_table_capacitance_rise[0][k]))
        minIndex1 = numpy.argmin(distance)
        distance[minIndex1] = 999999
        minIndex2 = numpy.argmin(distance)
        # print(minIndex1)
        distance.clear()
        xvals = [timing_table_capacitance_rise[0][minIndex1], timing_table_capacitance_rise[0][minIndex2]]
        yvals = [timing_table_transition_rise[0][findMiddle(timing_table_transition_rise[0])]]
        # print(yvals)
        fvals = [timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][minIndex1], timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][minIndex2]]

        # print(xvals)
        # print(yvals)
        # print(fvals)
        # print(fvals)
        print(interp2d(xvals, yvals, fvals))

        # cell_delay.append(cell_delay_things(CellwithLoad[i].cell,interp2d(xvals,findMiddle(timing_table_transition_rise[0]), fvals)))



