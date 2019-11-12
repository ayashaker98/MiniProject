#!/usr/bin/env python
# coding: utf-8

# In[171]:


from liberty.types import Group, select_cell, select_pin, select_timing_table
from liberty.parser import parse_liberty
from collections import namedtuple
import array
import numpy
from scipy.spatial.distance import euclidean
from scipy import arange, array, exp
from scipy.interpolate import interp1d, interp2d, InterpolatedUnivariateSpline
liberty_file = "/Users/noha/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())
fname = "/Users/noha/Desktop/rca4.rtlnopwr.v"

import re

maxFanoutin = 3 #input("Please Enter the maximum Fanout for any cell : ") 
maxDelay = 0.23 #input("Please Enter the maximum Delay for any cell : ") 
option = 3 #input("Please Enter the number of the method you would like to use: 1)Sizing Up.    2)Cloning.    3)Adding Buffers ") 
maxFanout = int(maxFanoutin)

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

#allpins has the information of every pin
for f in range(len(allpins)):
    cellgp = library.get_group('cell', allpins[f].cell.split('_')[0])
    pingp = cellgp.get_group('pin', allpins[f].pin)
    direct = pingp.__getitem__('direction')
    if (direct == 'input'):
        allpins[f] = allpins[f]._replace(direction = "input")
    else:
        allpins[f] = allpins[f]._replace(direction = "output")


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
        Fanout.append(FanoutCouple(allpins[l].cell,Fanoutnumber,"none")) #calculating the fanout of each cell before enhancement
        Fanoutnumber = -1

#checking wether or not the Fanout of each cell violates the maximum fanout
for g in range(len(Fanout)):
    if (Fanout[g].fanout > maxFanout):
        Fanout[g] = Fanout[g]._replace(status = "Violating")
    else:
        Fanout[g] = Fanout[g]._replace(status = "Not Violating")


finalpins = namedtuple('finalpins', ['cellname', 'capacitance','pinwithmaxcap'])
CapPins = []
MaxCapPins = []

#finding the capcitances of the input pins of each cell
for i in range(len(Fanout)):
    for j in range(len(allpins)):
        if(Fanout[i].cell == allpins[j].cell):
            cgroup = library.get_group('cell', allpins[j].cell.split('_')[0])
            pingroup = cgroup.get_group('pin', allpins[j].pin)
            cap = pingroup.__getitem__('capacitance')
            CapPins.append(finalpins(allpins[j].cell, cap, allpins[j].pin))
    MaxCapPins.append(finalpins(CapPins[0].cellname, max(CapPins, key=lambda k: k.capacitance).capacitance, max(CapPins, key=lambda k: k.capacitance).pinwithmaxcap))
   # CapPins.clear()

##finding the load capacitance of each fanout of each cell seperately
import collections
Fanout2 = collections.namedtuple('Fanoutpairs', 'gate1 gate2 outputpin1 inputpin2 pin1 pin2 inputcap')
Fanoutpairs = []
for i in range(len(allpins)):
    for j in range(len(allpins)):
        if ((allpins[i].parameter  == allpins[j].parameter)&(allpins[i].direction == 'output')&(allpins[j].direction =='input')):
            Fanoutpairs.append(Fanout2(allpins[i].cell,allpins[j].cell,allpins[i].parameter,allpins[j].parameter,allpins[i].pin,allpins[j].pin,0))

##finding the load capacitance of each fanout of each cell seperatelys
for i in range(len(Fanoutpairs)):
    for j in range(len(CapPins)):
        if((CapPins[j].cellname == Fanoutpairs[i].gate2) & (CapPins[j].pinwithmaxcap == Fanoutpairs[i].pin2)):
            Fanoutpairs[i]= Fanoutpairs[i]._replace(inputcap = CapPins[j].capacitance)

#summing all the capacitances of the fanouts of each cell --> CellwithLoad has the total load fanout of each cell
FanoutSum = 0
Fanout4 = namedtuple('Fanout4', ['cell', 'load', 'outputpin', 'relatedpin'])
CellwithLoad = []
for f in range(len(Fanout)):
    for r in range(len(Fanoutpairs)):
        if (Fanout[f].cell == Fanoutpairs[r].gate1):
            FanoutSum = FanoutSum + Fanoutpairs[r].inputcap
    CellwithLoad.append(Fanout4(Fanout[f].cell, FanoutSum, 'none', 'none'))
    FanoutSum = 0


#specifying the direction of each pin in CellwithLoad
for u in range(len(CellwithLoad)):
    for g in range(len(allpins)):
        if(CellwithLoad[u].cell == allpins[g].cell):
            if(allpins[g].direction == "output"):
                CellwithLoad[u] = CellwithLoad[u]._replace(outputpin = allpins[g].pin)
            elif((allpins[g].direction == "input") and (CellwithLoad[u].relatedpin == 'none')):
                CellwithLoad[u] = CellwithLoad[u]._replace(relatedpin = allpins[g].pin)


####################in this part we calculate the delay of each cell before any enahncement####################
rise_delay = -1
fall_delay = -1
minIndex1 = 0
minIndex2 = 0
o = 0
distance = []
distance2 = []
cell_delay_things = namedtuple('cell_delay_things', ['cell', 'delay', 'status','library'])
cell_delay_rise = []
cell_delay_fall = []
cell_delay = []
for i in range(len(CellwithLoad)):
    cellgroup = library.get_group('cell', CellwithLoad[i].cell.split('_')[0])
    pin = cellgroup.get_group('pin', CellwithLoad[i].outputpin)
    ##################timing tables for the rise delays################
    timing_table_rise = select_timing_table(pin, CellwithLoad[i].relatedpin, 'cell_rise')
    timing_table_transition_rise = timing_table_rise.get_array("index_1")
    timing_table_capacitance_rise = timing_table_rise.get_array("index_2")
    timing_table_values_rise = timing_table_rise.get_array("values")
    ##################timing tables for the fall delays#################
    timing_table_fall = select_timing_table(pin, CellwithLoad[i].relatedpin, 'cell_fall')
    timing_table_transition_fall = timing_table_fall.get_array("index_1")
    timing_table_capacitance_fall = timing_table_fall.get_array("index_2")
    timing_table_values_fall = timing_table_fall.get_array("values")
    ###################For calculating Rise delay#######################
    for j in range(len(timing_table_capacitance_rise[0])):  ###if the capacitance is in the timing table capacitances, we set it as the cell delay
        if (CellwithLoad[i].load == timing_table_capacitance_rise[0][j]):
            rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][j]

    if(CellwithLoad[i].load == 0):    ###if the output wire of the cell is an output to the medule we take the middle delay from the table
        rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])]
        fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])]
        cell_delay_rise.append(cell_delay_things(CellwithLoad[i].cell, timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])], 'none','none'))
        cell_delay_fall.append(cell_delay_things(CellwithLoad[i].cell, timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])], 'none','none'))

    if(rise_delay == -1): #####if the capacitance was not found in the table, then we must do interpolation or extrapolation
        xvals = timing_table_capacitance_rise[0]
        yvals = [timing_table_transition_rise[0][2]]
        fvals = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])]
        interp_rise = InterpolatedUnivariateSpline(xvals, fvals)
        cell_delay_rise.append(cell_delay_things(CellwithLoad[i].cell,interp_rise(CellwithLoad[i].load), 'none','none'))

    ###################For calculating Fall delay#######################
    for j in range(len(timing_table_capacitance_fall[0])):
        if (CellwithLoad[i].load == timing_table_capacitance_fall[0][j]):
            fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][j]

    if(fall_delay == -1):#####if the capacitance was not found in the table, then we must do interpolation or extrapolation
        xvals_fall = timing_table_capacitance_fall[0]
        yvals_fall = [timing_table_transition_fall[0][2]]
        fvals_fall = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])]
        interp_fall = InterpolatedUnivariateSpline(xvals_fall, fvals_fall)
        cell_delay_fall.append(cell_delay_things(CellwithLoad[i].cell,interp_fall(CellwithLoad[i].load), 'none','none'))
    rise_delay = -1
    fall_delay = -1


for q in range(len(cell_delay_fall)):  ##adding the maximum between rise and fall delays to cell_delay which has each cell with its delay
    maxFallRise = max(cell_delay_fall[q].delay, cell_delay_rise[q].delay)
    cell_delay.append(cell_delay_things(cell_delay_fall[q].cell, maxFallRise, 'none','none'))

for g in range(len(cell_delay)): ##checking if the delay is violating the maximum delay or not
    if (cell_delay[g].delay > float(maxDelay)):
        cell_delay[g] = cell_delay[g]._replace(status = "Violating")
    else:
        cell_delay[g] = cell_delay[g]._replace(status = "Not Violating")

with open('Bonus/original_delay.txt', 'w') as f:
    for item in range(len(cell_delay)):
        f.write(str(cell_delay[item].cell)+"     "+str(cell_delay[item].delay)+ " " + "\n")
################################################################################################################

numoftypes = namedtuple('numoftypes', ['cell', 'num'])
typeslist = []
numcelltype = 0
for k in range(len(CellwithLoad)):
    cellfixed = CellwithLoad[k].cell.split('_')[0]
    for z in range(len(CellwithLoad)):
        if(cellfixed == CellwithLoad[z].cell.split('_')[0]):
            numcelltype = numcelltype + 1
    typeslist.append(numoftypes(cellfixed, numcelltype))
    numcelltype = 0
typeslist = list(dict.fromkeys(typeslist))
# for p in range(len(typeslist)):
#     print(typeslist[p])

with open('Bonus/original_Number_of_cells.txt', 'w') as f:
    for item in range(len(typeslist)):
        f.write(str(typeslist[item].cell)+"     "+str(typeslist[item].num)+ " " + "\n")

#cell_delay has the final delay of each cell with a space saying if its violated or not
#Fanout has the final fanout of each cell with a space saying if its violated or not

import math as m
Bufferadded = namedtuple('Bufferadded', ['BufferType', 'Buffername','inputt','output'])
Buf = []
KeepTrack = namedtuple('KeepTrack', ['output', 'somethingelse','r','Stage','inputt','Fanout','done'])
keep = []
maxfanout = 2

def sizeup():
    f = open("/Users/noha/Desktop/osu035.lib", "r")
    librarycells = namedtuple('librarycells', ['cell', 'area'])
    libcell =[]
    allcells = []
    count = 0
    for x in f:
        if("cell (" in x):
            allcells.append(x.split('cell (')[1])
    allcells2 = namedtuple('allcells2', ['cell', 'area'])
    cells = []
    alllibcells = []
    for i in range(len(allcells)):
        alllibcells.append(allcells2(allcells[i].split(')')[0],0))
        cells.append(allcells[i].split(')')[0])  
    for i in range(len(alllibcells)):
        n = library.get_group('cell', cells[i])
        alllibcells[i] = alllibcells[i]._replace(area = n.__getitem__('area'))
    for i in range(len(cell_delay)):
        cell_delay[i] =  cell_delay[i]._replace(library =cell_delay[i].cell.split('_')[0])
    replaceit = namedtuple('replaceit', ['origcell','origgate','num','newcell','newarea','library'])
    repl = []
    for r in range(len(cell_delay)):
        repl.append(replaceit(cell_delay[r].cell,cell_delay[r].library[0:len(cell_delay[r].library)-1],cell_delay[r].library[len(cell_delay[r].library)-1],'none',0,cell_delay[r].cell.split('_')[0]))
    for j in range(len(repl)):
        if(cell_delay[j].status == 'Violating'):
            for x in range(len(alllibcells)):
                if((repl[j].origgate == alllibcells[x].cell[0:len(alllibcells[x].cell)-1])&(repl[j].num<alllibcells[x].cell[len(alllibcells[x].cell)-1])):
                    repl[j] = repl[j]._replace(newcell=alllibcells[x].cell,newarea =alllibcells[x].area)
                    
    for u in range(len(functions)):
        for h in range(len(repl)):
            if((repl[h].library in functions[u])&(repl[h].newcell!='none')):
                functions[u] = functions[u].replace(repl[h].library,repl[h].newcell)

    numoftypessu = namedtuple('numoftypes', ['cell', 'num'])
    typeslistsu = []
    numcelltype = 0
    for k in range(len(CellwithLoad)):
        cellfixed = CellwithLoad[k].cell.split('_')[0]
        for z in range(len(CellwithLoad)):
            if(cellfixed == CellwithLoad[z].cell.split('_')[0]):
                numcelltype = numcelltype + 1
        typeslistsu.append(numoftypessu(cellfixed, numcelltype))
        numcelltype = 0
    typeslistsu = list(dict.fromkeys(typeslistsu))

    with open('Bonus/Number_of_cells_after_sizing_up.txt', 'w') as f:
        for item2 in range(len(typeslistsu)):
            f.write(str(typeslistsu[item2].cell)+"     "+str(typeslistsu[item2].num)+ " " + "\n")

    return repl


# In[579]:

def addBuffers():
    maxim = maxFanout
    count = 0
    Stage = 0
    numBuf = 0
    n = 0
    r = []
    for i in range(len(Fanout)):
        maxim=maxFanout
        numBuf = m.ceil(Fanout[i].fanout/maxFanout)
        nmbf = numBuf
        flag = isinstance(numBuf, float)
        if(Fanout[i].status == 'Violating'): 
            if(numBuf<=maxFanout):
                Stage = 1
                while(nmbf!=0):
                    for p in range(len(allpins)):########Find input for first stage
                                    if((allpins[p].cell == Fanout[i].cell)&(allpins[p].direction=='output')):
                                        Buf.append(Bufferadded('BUFX4 ',"BUFX4_"+ str(count)+"_added " ,allpins[p].parameter,str(count) ))

                                        keep.append(KeepTrack(str(count),"none","none",Stage,0,0,'none'))
                                        count = count +1
                                        nmbf=nmbf-1
            else:
                n = Fanout[i].fanout/maximum
                while(n>1):
                    r.append(int(n))
                    n = (n/maxFanout)
                    r.reverse()
                    for z in range(len(r)):########Stage########
                        if(z==0):########Stage 1########
                                Stage = 1
                                for w in range(len(allpins)):########Find input for first stage  
                                     for k in range(r[z]):####Add buffers according to r[i] in first stage########## 
                                        if((allpins[w].cell == Fanout[i].cell)&(allpins[w].direction=='output')):
                                            Buf.append(Bufferadded('BUFX4 ',"BUFX4_"+ str(count)+"_added " ,allpins[w].parameter,str(count)))
                                            keep.append(KeepTrack(str(count),z,r[0],Stage,allpins[w].parameter,Fanout[i].cell,'none'))
                                            count = count + 1

                                            break
                        else:######MultiStage
                            Stage = Stage + 1
                            for y in range(r[z]):
                                Buf.append(Bufferadded('BUFX4 ',"BUFX4_"+ str(count)+"_added " ,Buf[len(Buf)-1].output,str(count)))
                                keep.append(KeepTrack(str(count+1),z,r[0],Stage,str(count),Fanout[i].cell,'none'))
                                count = count+1                      
        n = Stage
        jj = numBuf
        q = len(r)
        counter=0
        e = 0
        if(len(keep)<numBuf):
            e = 1
        elif(len(keep)==maxFanout):
            e=0
        elif(maxFanout%2==0):
            e = maxFanout
        elif(len(Buf)%2==0):
            e=1

        for f in range(len(keep)):
            if((keep[f-e].r != 'none')&(keep[f-e].Stage == 1)&(keep[f].Stage==2)&(keep[f].done=='none')):
                while(counter!=n):
                    while(jj!=0):
                        Buf[f] = Buf[f]._replace(inputt = Buf[f-e].output)
                        keep[f] = keep[f]._replace(done = "done")

                        jj=jj-1
                    jj=e 
                    counter = counter+1
            elif ((keep[f-e].r!='none')&(keep[f-e].Stage==2)&(keep[f].Stage==2)&(keep[f-e].done=='done')):
                while(counter!=n):
                    while(jj!=0):
                        Buf[f] = Buf[f]._replace(inputt = Buf[f-e].inputt)
                        keep[f] = keep[f]._replace(done = "done")
                        jj=jj-1
                    jj=e 
                    counter = counter+1
            counter = 0        
        r.clear()  
    corrpin = namedtuple('corrpin', ['gateoutput', 'finalout','maxfanout'])
    cor = []
    Finalpin = namedtuple('Finalpin', ['output', 'input'])
    fin = []
    fin2 = []
    flag = 0
    for d in range(len(Buf)):
        fin.append(Buf[d].inputt)
        fin2.append(Buf[d].output)
    for x in range(len(Buf)):  
        n = str(Buf[x].output)
        g = str(Buf[x].inputt)
        if((n in fin) == False)&((g in fin2)==False):
            cor.append(corrpin(Buf[x].inputt,Buf[x].output,maxFanout))
        for y in range(len(Buf)):
            if((Buf[x].output == Buf[y].inputt)):
                cor.append(corrpin(Buf[x].inputt,Buf[y].output,maxFanout))
    counting = 0
    for y in range(len(functions)):
        counting = 0
        for u in range(len(cor)):
            if(".Y("+cor[u].gateoutput+")" in functions[y]):
                break
            elif((cor[u].gateoutput in functions[y])&(cor[u].maxfanout>0)):
                functions[y] = functions[y].replace(cor[u].gateoutput,str(cor[u].finalout))
                l = cor[u].maxfanout-1
                cor[u] = cor[u]._replace(maxfanout =l)

    numoftypesbuf = namedtuple('numoftypes', ['cell', 'num'])
    typeslistbuf = []
    numcelltype = 0
    numcelltypebuf = 0
    flagfound = False
    for k in range(len(CellwithLoad)):
        cellfixed = CellwithLoad[k].cell.split('_')[0]
        for z in range(len(CellwithLoad)):
            if(cellfixed == CellwithLoad[z].cell.split('_')[0]):
                numcelltype = numcelltype + 1
        typeslistbuf.append(numoftypesbuf(cellfixed, numcelltype))
        numcelltype = 0
    typeslistbuf = list(dict.fromkeys(typeslistbuf))

    for f in range(len(Buf)):
        numcelltypebuf = numcelltypebuf + 1
    
    for l in range(len(typeslistbuf)):
        if(typeslistbuf[l].cell == Buf[0].BufferType):
            typeslistbuf[l] = typeslistbuf[l]._repalce(num = num + 1)
            flagfound = True
    if(flagfound == False):
        typeslistbuf.append(numoftypesbuf(Buf[0].BufferType, numcelltypebuf))

    with open('Bonus/Number_of_cells_after_adding_buffers.txt', 'w') as f:
        for item3 in range(len(typeslistbuf)):
            f.write(str(typeslistbuf[item3].cell)+"     "+str(typeslistbuf[item3].num)+ " " + "\n")
    flagfound = False
    # print(Buf)
    return Buf

###########################Cloning#########################################
def Cloning(functions, dictList, allpins, Fanout):
    newNetlist = []
    paramstuple = namedtuple('paramstuple', ['cell','pp', 'pin', 'param', 'direction'])
    clonesstr = ""
    allclonedcells = []
    originalclone = []
    inoutstr = ""
    paramsCounter = 0
    startp = '.'
    endp = ','
    paramslist = []

    for r in range(len(Fanout)):
        if(Fanout[r].status == "Violating"):
            numClones = Fanout[r].fanout/int(maxFanout)
            if(numClones.is_integer() == False):   ###specifying how many clones I need
                numClones = int(numClones) + 1 
            else:
                numClones = int(numClones)
            for x in range(numClones-1):
                for k in range(len(dictList)):
                    if(Fanout[r].cell == dictList[k][0]):
                        inoutstr = inoutstr + "("
                        for y in range(len(dictList[k][1])):
                            inoutstr = inoutstr + dictList[k][1][y]
                        inoutstr = inoutstr + ")"
                        clonesstr = clonesstr + Fanout[r].cell.split('_')[0]+" "+Fanout[r].cell+"_cloned"+str(r)+ " "+ inoutstr
                        originalclone.append(clonesstr)   ##cloning the violating cells
                        inoutstr = ""
                libcell = clonesstr.split(' ')[0]
                cellinstance = clonesstr.split(' ')[1]
                params = clonesstr.split(' ')[2]
                for o in range(len(params)):
                    if(params[o] == '.'):
                        paramsCounter = paramsCounter + 1
                for p in range(paramsCounter):
                    paramslist.append(paramstuple(cellinstance,params[params.find(startp)+len(startp):params.find(endp)], 'none', 'none', 'none'))
                    params = params.replace(params[params.find(startp)+len(startp):params.find(endp)], '')
                    params = params.replace(startp, '', 1)
                    params = params.replace(endp, '', 1)
                    
                clonesstr = ""
                paramsCounter = 0
    ##analyzing the cloned cells and specifying their pins and parameters
    ss = '('
    ee = ')'
    for i in range(len(paramslist)):
        pp2 = paramslist[i].pp
        paramslist[i] = paramslist[i]._replace(pin = pp2.split('(')[0])
        paramslist[i] = paramslist[i]._replace(param = pp2[pp2.find(ss)+len(ss):pp2.find(ee)])
        cellgp2 = library.get_group('cell', paramslist[i].cell.split('_')[0])
        pingp2 = cellgp2.get_group('pin', paramslist[i].pin)
        direct2 = pingp2.__getitem__('direction')
        paramslist[i] = paramslist[i]._replace(direction = direct2)
        if(paramslist[i].direction == 'output'):
            paramslist[i] = paramslist[i]._replace(param = pp2[pp2.find(ss)+len(ss):pp2.find(ee)]+paramslist[i].param[-1:])

    
    allpins_clones = []
    allpinsforclones = []
    otherallpins = allpins
    pins_clones = namedtuple('pins', ['cell','pin', 'parameter', 'direction', 'cloneflag'])
    for k in range(len(allpins)):
        allpins_clones.append(pins_clones(allpins[k].cell, allpins[k].pin, allpins[k].parameter, allpins[k].direction, False))

    for t in range(len(paramslist)):
        allpins_clones.append(pins_clones(paramslist[t].cell, paramslist[t].pin, paramslist[t].param, paramslist[t].direction, True))
        otherallpins.append(pins_clones(paramslist[t].cell, paramslist[t].pin, paramslist[t].param, paramslist[t].direction, True))

    ##changing the input of a cell to be the output of one of the clones
    for j in range(len(allpins_clones)):  
        if(allpins_clones[j].cloneflag == True):
            for l in range(len(allpins)):
                if((allpins_clones[j].cell.split('_cloned')[0] == allpins[l].cell) and (allpins_clones[j].direction == 'output') and (allpins_clones[j].parameter[0:-1] == allpins[l].parameter)):
                    for h in range(maxFanout):
                        if(allpins[l].direction == 'output'):
                            for v in range(len(allpins)):
                                if((allpins[v].parameter == allpins[l].parameter) and (allpins[v].direction == 'input')):
                                    otherallpins[v] = allpins[v]._replace(parameter = allpins_clones[j].parameter)
                                    allpinsforclones.append(allpins[v]._replace(parameter = allpins_clones[j].parameter))
                                    break

    allpins_final_clones = []
    pins_clones_final = namedtuple('pins', ['cell','pin', 'parameter', 'direction'])

    for u in range(len(paramslist)):
        allpins_final_clones.append(pins_clones_final(paramslist[u].cell, paramslist[u].pin, paramslist[u].param, paramslist[u].direction))

    ###recreating the netlist
    for z in range(len(functions)):
        for c in range(len(allpinsforclones)):
            if(allpinsforclones[c].cell in functions[z]):
                functions[z] = functions[z].replace(allpinsforclones[c].parameter[0:-1], allpinsforclones[c].parameter)

    for r in range(len(originalclone)):
        for c in range(len(paramslist)):
            if(paramslist[c].cell in originalclone[r]):
                if (paramslist[c].direction == 'output'):
                    originalclone[r] = originalclone[r].replace(paramslist[c].param[0:-1], paramslist[c].param)
                else:
                    originalclone[r] = originalclone[r].replace(paramslist[c].param, paramslist[c].param)


    for g in range(len(originalclone)):
        functions.append(originalclone[g])
    
    #####finding the number of cells of each type for the bonus###########################
    numoftypescl = namedtuple('numoftypes', ['cell', 'num'])
    typeslistcl = []
    numcelltype = 0
    for k in range(len(functions)):
        cellfixed = functions[k].split(' ')[0]
        for z in range(len(functions)):
            if(cellfixed == functions[z].split(' ')[0]):
                numcelltype = numcelltype + 1
        typeslistcl.append(numoftypescl(cellfixed, numcelltype))
        numcelltype = 0
    typeslistcl = list(dict.fromkeys(typeslistcl))

    with open('Bonus/Number_of_cells_after_cloning.txt', 'w') as f:
        for item4 in range(len(typeslistcl)):
            f.write(str(typeslistcl[item4].cell)+"     "+str(typeslistcl[item4].num)+ " " + "\n")

    #################recalculating the delay after enhancement for the bonus##############
    Fanoutcl = collections.namedtuple('Fanoutpairs', 'gate1 gate2 outputpin1 inputpin2 pin1 pin2 inputcap')
    Fanoutpairscl = []
    for i in range(len(otherallpins)):
        for j in range(len(otherallpins)):
            if ((otherallpins[i].parameter  == otherallpins[j].parameter)&(otherallpins[i].direction == 'output')&(otherallpins[j].direction =='input')):
                Fanoutpairscl.append(Fanoutcl(otherallpins[i].cell,otherallpins[j].cell,otherallpins[i].parameter,otherallpins[j].parameter,otherallpins[i].pin,otherallpins[j].pin,0))

    finalpinscl = namedtuple('finalpins', ['cellname', 'capacitance','pinwithmaxcap'])
    CapPinscl = []

    for i in range(len(functions)):
        for j in range(len(otherallpins)):
            if(functions[i].split(' ')[1] == otherallpins[j].cell):
                cgroup = library.get_group('cell', otherallpins[j].cell.split('_')[0])
                pingroup = cgroup.get_group('pin', otherallpins[j].pin)
                cap = pingroup.__getitem__('capacitance')
                CapPinscl.append(finalpinscl(otherallpins[j].cell, cap, otherallpins[j].pin))

    for i in range(len(Fanoutpairscl)):
        for j in range(len(CapPinscl)):
            if((CapPinscl[j].cellname == Fanoutpairscl[i].gate2) & (CapPinscl[j].pinwithmaxcap == Fanoutpairscl[i].pin2)):
                Fanoutpairscl[i]= Fanoutpairscl[i]._replace(inputcap = CapPinscl[j].capacitance)

    FanoutSumcl = 0
    Fanout4cl = namedtuple('Fanout4', ['cell', 'load', 'outputpin', 'relatedpin'])
    CellwithLoadcl = []
    for f in range(len(functions)):
        for r in range(len(Fanoutpairscl)):
            if (functions[f].split(' ')[1] == Fanoutpairscl[r].gate1):
                FanoutSumcl = FanoutSumcl + Fanoutpairscl[r].inputcap
        CellwithLoadcl.append(Fanout4cl(functions[f].split(' ')[1], FanoutSumcl, 'none', 'none'))
        FanoutSumcl = 0

    for u in range(len(CellwithLoadcl)):
        for g in range(len(otherallpins)):
            if(CellwithLoadcl[u].cell == otherallpins[g].cell):
                if(otherallpins[g].direction == "output"):
                    CellwithLoadcl[u] = CellwithLoadcl[u]._replace(outputpin = otherallpins[g].pin)
                elif((otherallpins[g].direction == "input") and (CellwithLoadcl[u].relatedpin == 'none')):
                    CellwithLoadcl[u] = CellwithLoadcl[u]._replace(relatedpin = otherallpins[g].pin)

    rise_delay = -1
    fall_delay = -1
    minIndex1 = 0
    minIndex2 = 0
    o = 0
    distance = []
    distance2 = []
    cell_delay_things = namedtuple('cell_delay_things', ['cell', 'delay', 'status','library'])
    cell_delay_rise = []
    cell_delay_fall = []
    cell_delay = []
    for i in range(len(CellwithLoadcl)):
        cellgroup = library.get_group('cell', CellwithLoadcl[i].cell.split('_')[0])
        pin = cellgroup.get_group('pin', CellwithLoadcl[i].outputpin)
        ##################timing tables for the rise delays################
        timing_table_rise = select_timing_table(pin, CellwithLoadcl[i].relatedpin, 'cell_rise')
        timing_table_transition_rise = timing_table_rise.get_array("index_1")
        timing_table_capacitance_rise = timing_table_rise.get_array("index_2")
        timing_table_values_rise = timing_table_rise.get_array("values")
        ##################timing tables for the fall delays#################
        timing_table_fall = select_timing_table(pin, CellwithLoadcl[i].relatedpin, 'cell_fall')
        timing_table_transition_fall = timing_table_fall.get_array("index_1")
        timing_table_capacitance_fall = timing_table_fall.get_array("index_2")
        timing_table_values_fall = timing_table_fall.get_array("values")
        ###################For calculating Rise delay#######################
        for j in range(len(timing_table_capacitance_rise[0])):
            if (CellwithLoadcl[i].load == timing_table_capacitance_rise[0][j]):
                rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][j]
                print("Rise Delay was set successfuly")

        if(CellwithLoadcl[i].load == 0):     ###if the output wire of the cell is an output to the medule we take the middle delay from the table
            rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])]
            fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])]
            cell_delay_rise.append(cell_delay_things(CellwithLoadcl[i].cell, timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])], 'none','none'))
            cell_delay_fall.append(cell_delay_things(CellwithLoadcl[i].cell, timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])], 'none','none'))

        if(rise_delay == -1):      #####if the capacitance was not found in the table, then we must do interpolation or extrapolation
            xvals = timing_table_capacitance_rise[0]
            yvals = [timing_table_transition_rise[0][2]]
            fvals = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])]
            interp_rise = InterpolatedUnivariateSpline(xvals, fvals)
            cell_delay_rise.append(cell_delay_things(CellwithLoadcl[i].cell,interp_rise(CellwithLoadcl[i].load), 'none','none'))

        ###################For calculating Fall delay#######################
        for j in range(len(timing_table_capacitance_fall[0])):
            if (CellwithLoadcl[i].load == timing_table_capacitance_fall[0][j]):
                fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][j]
                print("Fall Delay was set successfuly")

        if(fall_delay == -1):  #####if the capacitance was not found in the table, then we must do interpolation or extrapolation
            xvals_fall = timing_table_capacitance_fall[0]
            yvals_fall = [timing_table_transition_fall[0][2]]
            fvals_fall = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])]
            interp_fall = InterpolatedUnivariateSpline(xvals_fall, fvals_fall)
            cell_delay_fall.append(cell_delay_things(CellwithLoadcl[i].cell,interp_fall(CellwithLoadcl[i].load), 'none','none'))
        rise_delay = -1
        fall_delay = -1


    for q in range(len(cell_delay_fall)):
        maxFallRise = max(cell_delay_fall[q].delay, cell_delay_rise[q].delay)
        cell_delay.append(cell_delay_things(cell_delay_fall[q].cell, maxFallRise, 'none','none'))

    for g in range(len(cell_delay)):
        if (cell_delay[g].delay > float(maxDelay)):
            cell_delay[g] = cell_delay[g]._replace(status = "Violating")
        else:
            cell_delay[g] = cell_delay[g]._replace(status = "Not Violating")

        ###############printing the recalculated delay to a file###############
        with open('Bonus/Cell_delay_after_cloning.txt', 'w') as f:
            for item5 in range(len(cell_delay)):
                f.write(str(cell_delay[item5].cell)+"     "+str(cell_delay[item5].delay)+  "\n")

    return functions
##################################Cloning######################################


##########PICK##################
if(int(option) ==1):
    n = sizeup()    ##do sizing up

    #######Recalculating the delay after sizing up for the bonus########
    allpins2 = allpins.copy()
    rrrr = namedtuple('rrrr', 'originalgate newcell newgate library inputcap')
    rr = []
    for u in range(len(n)):
        rr.append(rrrr(n[u].origcell,n[u].newcell,'',n[u].library,0))
    p=''
    for b in range(len(n)):
        p = rr[b].originalgate
        p = p.replace(rr[b].library,rr[b].newcell)
        rr[b]=rr[b]._replace(newgate = p)

    for t in range(len(allpins2)):
        for y in range(len(rr)):
            if((rr[y].originalgate == allpins2[t].cell)&(rr[y].newcell!='none')):
                allpins2[t]=allpins2[t]._replace(cell = rr[y].newgate)

    part2 = namedtuple('part1', 'gate1 gate2 outputpin1 inputpin2 pin1 pin2 inputcap')
    Fanoutpairs2 = []
    for i in range(len(allpins2)):
        for j in range(len(allpins2)):
            if ((allpins2[i].parameter  == allpins2[j].parameter)&(allpins2[i].direction == 'output')&(allpins2[j].direction =='input')):
                Fanoutpairs2.append(part2(allpins2[i].cell,allpins2[j].cell,allpins2[i].parameter,allpins2[j].parameter,allpins2[i].pin,allpins2[j].pin,0))

    finalpinscl = namedtuple('finalpins', ['cellname', 'capacitance','pinwithmaxcap'])
    CapPinssu = []
    
    for i in range(len(functions)):
        for j in range(len(allpins2)):
            if(functions[i].split(' ')[1] == allpins2[j].cell):
                cgroup = library.get_group('cell', allpins2[j].cell.split('_')[0])
                pingroup = cgroup.get_group('pin', allpins2[j].pin)
                cap = pingroup.__getitem__('capacitance')
                CapPinssu.append(finalpinscl(allpins2[j].cell, cap, allpins2[j].pin))
      
    for i in range(len(Fanoutpairs2)):
        for j in range(len(CapPinssu)):
            if((CapPinssu[j].cellname == Fanoutpairs2[i].gate2) & (CapPinssu[j].pinwithmaxcap == Fanoutpairs2[i].pin2)):
                Fanoutpairs2[i]= Fanoutpairs2[i]._replace(inputcap = CapPinssu[j].capacitance)
    
    FanoutSumsu = 0
    Fanout4su = namedtuple('Fanout4', ['cell', 'load', 'outputpin', 'relatedpin'])
    CellwithLoadsu = []
    for f in range(len(functions)):
        for r in range(len(Fanoutpairs2)):
            if (functions[f].split(' ')[1] == Fanoutpairs2[r].gate1):
                FanoutSumsu = FanoutSumsu + Fanoutpairs2[r].inputcap
        CellwithLoadsu.append(Fanout4su(functions[f].split(' ')[1], FanoutSumsu, 'none', 'none'))
        FanoutSumsu = 0

    for u in range(len(CellwithLoadsu)):
        for g in range(len(allpins2)):
            if(CellwithLoadsu[u].cell == allpins2[g].cell):
                if(allpins2[g].direction == "output"):
                    CellwithLoadsu[u] = CellwithLoadsu[u]._replace(outputpin = allpins2[g].pin)
                elif((allpins2[g].direction == "input") and (CellwithLoadsu[u].relatedpin == 'none')):
                    CellwithLoadsu[u] = CellwithLoadsu[u]._replace(relatedpin = allpins2[g].pin)
        
    rise_delay = -1
    fall_delay = -1
    minIndex1 = 0
    minIndex2 = 0
    o = 0
    distance = []
    distance2 = []
    cell_delay_things = namedtuple('cell_delay_things', ['cell', 'delay', 'status','library'])
    cell_delay_rise = []
    cell_delay_fall = []
    cell_delay = []
    for i in range(len(CellwithLoadsu)):
        cellgroup = library.get_group('cell', CellwithLoadsu[i].cell.split('_')[0])
        pin = cellgroup.get_group('pin', CellwithLoadsu[i].outputpin)
        ##################timing tables for the rise delays################
        timing_table_rise = select_timing_table(pin, CellwithLoadsu[i].relatedpin, 'cell_rise')
        timing_table_transition_rise = timing_table_rise.get_array("index_1")
        timing_table_capacitance_rise = timing_table_rise.get_array("index_2")
        timing_table_values_rise = timing_table_rise.get_array("values")
        ##################timing tables for the fall delays#################
        timing_table_fall = select_timing_table(pin, CellwithLoadsu[i].relatedpin, 'cell_fall')
        timing_table_transition_fall = timing_table_fall.get_array("index_1")
        timing_table_capacitance_fall = timing_table_fall.get_array("index_2")
        timing_table_values_fall = timing_table_fall.get_array("values")
        ###################For calculating Rise delay#######################
        for j in range(len(timing_table_capacitance_rise[0])):
            if (CellwithLoadsu[i].load == timing_table_capacitance_rise[0][j]):
                rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][j]
                print("Rise Delay was set successfuly")

        if(CellwithLoadsu[i].load == 0):   ##if the output wire of the cell is an output to the medule we take the middle delay from the table
            rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])]
            fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])]
            cell_delay_rise.append(cell_delay_things(CellwithLoadsu[i].cell, timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])], 'none','none'))
            cell_delay_fall.append(cell_delay_things(CellwithLoadsu[i].cell, timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])], 'none','none'))
 
        if(rise_delay == -1):     #####if the capacitance was not found in the table, then we must do interpolation or extrapolation
            xvals = timing_table_capacitance_rise[0]
            yvals = [timing_table_transition_rise[0][2]]
            fvals = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])]
            interp_rise = InterpolatedUnivariateSpline(xvals, fvals)
            cell_delay_rise.append(cell_delay_things(CellwithLoadsu[i].cell,interp_rise(CellwithLoadsu[i].load), 'none','none'))

        ###################For calculating Fall delay#######################
        for j in range(len(timing_table_capacitance_fall[0])):
            if (CellwithLoadsu[i].load == timing_table_capacitance_fall[0][j]):
                fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][j]
                print("Fall Delay was set successfuly")

        if(fall_delay == -1):     #####if the capacitance was not found in the table, then we must do interpolation or extrapolation  
            xvals_fall = timing_table_capacitance_fall[0]
            yvals_fall = [timing_table_transition_fall[0][2]]
            fvals_fall = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])]
            interp_fall = InterpolatedUnivariateSpline(xvals_fall, fvals_fall)
            cell_delay_fall.append(cell_delay_things(CellwithLoadsu[i].cell,interp_fall(CellwithLoadsu[i].load), 'none','none'))
        rise_delay = -1
        fall_delay = -1


    for q in range(len(cell_delay_fall)):
        maxFallRise = max(cell_delay_fall[q].delay, cell_delay_rise[q].delay)
        cell_delay.append(cell_delay_things(cell_delay_fall[q].cell, maxFallRise, 'none','none'))

    for g in range(len(cell_delay)):
        if (cell_delay[g].delay > float(maxDelay)):
            cell_delay[g] = cell_delay[g]._replace(status = "Violating")
        else:
            cell_delay[g] = cell_delay[g]._replace(status = "Not Violating")

    #############printing the recalculated delay################
    with open('Bonus/Cell_delay_after_sizing_up.txt', 'w') as f:
        for item6 in range(len(cell_delay)):
            f.write(str(cell_delay[item6].cell)+"     "+str(cell_delay[item6].delay)+  "\n")

elif(int(option) == 2):
    cloned = Cloning(functions, dictList, allpins, Fanout)     ####performing cloning####
else:
    n = addBuffers()
    # print(n)
    Bufout = []
    for j in range(len(n)):
        Bufout.append(n[j].BufferType+Buf[j].Buffername+" ( .A("+n[j].inputt+"), "+".Y("+n[j].output+"));")
    functions = functions + Bufout
    dictList2 = []
    for code in functions:
        n2 = code.split(maxsplit=2)
        del n2[0]
        dictList2.append(n2)
        
    ##############################recalculating the delay after adding buffers for the bonus##################
    for i in range(len(dictList2)):
        for j in range(2):
            n2 =  re.sub(";", '', dictList2[i][j])
            dictList2[i][j] = n2


            
    for i in range(len(dictList2)):
        n2 = dictList2[i][1].split()
        dictList2[i][1] = n2
    for i in range(len(dictList2)):
        del dictList2[i][1][0]
    start = '.'
    end = '('
    end2 = ')'
    pins2 = namedtuple('pins', ['cell','pin', 'parameter', 'direction'])
    allpins3 = []

    for x in range(len(dictList2)):
        for y in range(len(dictList2[x])):
            for z in range(len(dictList2[x][y])):
                if (dictList2[x][y][z][dictList2[x][y][z].find(start)+len(start):dictList2[x][y][z].rfind(end)] != ""):
                    allpins3.append(pins(dictList2[x][0], dictList2[x][y][z][dictList2[x][y][z].find(start)+len(start):dictList2[x][y][z].rfind(end)], dictList2[x][y][z][dictList2[x][y][z].find(end)+len(end):dictList2[x][y][z].rfind(end2)], "none"))


    for f in range(len(allpins3)):
        cellgp = library.get_group('cell', allpins3[f].cell.split('_')[0])
        pingp = cellgp.get_group('pin', allpins3[f].pin)
        direct = pingp.__getitem__('direction')
        if (direct == 'input'):
            allpins3[f] = allpins3[f]._replace(direction = "input")
        else:
            allpins3[f] = allpins3[f]._replace(direction = "output")
    Fanout3 = namedtuple('Fanoutpairs', 'gate1 gate2 outputpin1 inputpin2 pin1 pin2 inputcap')
    Fanoutpairs3 = []
    for i in range(len(allpins3)):
        for j in range(len(allpins3)):
            if ((allpins3[i].parameter  == allpins3[j].parameter)&(allpins3[i].direction == 'output')&(allpins3[j].direction =='input')):
                Fanoutpairs3.append(Fanout3(allpins3[i].cell,allpins3[j].cell,allpins3[i].parameter,allpins3[j].parameter,allpins3[i].pin,allpins3[j].pin,0))

    
    finalpinsbf = namedtuple('finalpins', ['cellname', 'capacitance','pinwithmaxcap'])
    CapPinsbf = []
    
    for i in range(len(functions)):
        for j in range(len(allpins3)):
            if(functions[i].split(' ')[1] == allpins3[j].cell):
                cgroup = library.get_group('cell', allpins3[j].cell.split('_')[0])
                pingroup = cgroup.get_group('pin', allpins3[j].pin)
                cap = pingroup.__getitem__('capacitance')
                CapPinsbf.append(finalpinsbf(allpins3[j].cell, cap, allpins3[j].pin))
    
    
    for i in range(len(Fanoutpairs3)):
        for j in range(len(CapPinsbf)):
            if((CapPinsbf[j].cellname == Fanoutpairs3[i].gate2) & (CapPinsbf[j].pinwithmaxcap == Fanoutpairs3[i].pin2)):
                Fanoutpairs3[i]= Fanoutpairs3[i]._replace(inputcap = CapPinsbf[j].capacitance)

    FanoutSumbf = 0
    Fanout4bf = namedtuple('Fanout4', ['cell', 'load', 'outputpin', 'relatedpin'])
    CellwithLoadbf = []
    for f in range(len(functions)):
        for r in range(len(Fanoutpairs3)):
            if (functions[f].split(' ')[1] == Fanoutpairs3[r].gate1):
                FanoutSumbf = FanoutSumbf + Fanoutpairs3[r].inputcap
        CellwithLoadbf.append(Fanout4bf(functions[f].split(' ')[1], FanoutSumbf, 'none', 'none'))
        FanoutSumbf = 0

    for u in range(len(CellwithLoadbf)):
        for g in range(len(allpins3)):
            if(CellwithLoadbf[u].cell == allpins3[g].cell):
                if(allpins3[g].direction == "output"):
                    CellwithLoadbf[u] = CellwithLoadbf[u]._replace(outputpin = allpins3[g].pin)
                elif((allpins3[g].direction == "input") and (CellwithLoadbf[u].relatedpin == 'none')):
                    CellwithLoadbf[u] = CellwithLoadbf[u]._replace(relatedpin = allpins3[g].pin)
        
    rise_delay = -1
    fall_delay = -1
    minIndex1 = 0
    minIndex2 = 0
    o = 0
    distance = []
    distance2 = []
    cell_delay_things = namedtuple('cell_delay_things', ['cell', 'delay', 'status','library'])
    cell_delay_rise = []
    cell_delay_fall = []
    cell_delay = []
    for i in range(len(CellwithLoadbf)):
        cellgroup = library.get_group('cell', CellwithLoadbf[i].cell.split('_')[0])
        pin = cellgroup.get_group('pin', CellwithLoadbf[i].outputpin)
        ##################timing tables for the rise delays################
        timing_table_rise = select_timing_table(pin, CellwithLoadbf[i].relatedpin, 'cell_rise')
        timing_table_transition_rise = timing_table_rise.get_array("index_1")
        timing_table_capacitance_rise = timing_table_rise.get_array("index_2")
        timing_table_values_rise = timing_table_rise.get_array("values")
        ##################timing tables for the fall delays#################
        timing_table_fall = select_timing_table(pin, CellwithLoadbf[i].relatedpin, 'cell_fall')
        timing_table_transition_fall = timing_table_fall.get_array("index_1")
        timing_table_capacitance_fall = timing_table_fall.get_array("index_2")
        timing_table_values_fall = timing_table_fall.get_array("values")
        ###################For calculating Rise delay#######################
        for j in range(len(timing_table_capacitance_rise[0])):     ###if the capacitance is in the timing table capacitances, we set it as the cell delay
            if (CellwithLoadbf[i].load == timing_table_capacitance_rise[0][j]):
                rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][j]

        if(CellwithLoadbf[i].load == 0):   ###if the output wire of the cell is an output to the medule we take the middle delay from the table
            rise_delay = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])]
            fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])]
            cell_delay_rise.append(cell_delay_things(CellwithLoadbf[i].cell, timing_table_values_rise[findMiddle(timing_table_transition_rise[0])][findMiddle(timing_table_capacitance_rise[0])], 'none','none'))
            cell_delay_fall.append(cell_delay_things(CellwithLoadbf[i].cell, timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][findMiddle(timing_table_capacitance_fall[0])], 'none','none'))

        if(rise_delay == -1):     #####if the capacitance was not found in the table, then we must do interpolation or extrapolation  
            xvals = timing_table_capacitance_rise[0]
            yvals = [timing_table_transition_rise[0][2]]
            fvals = timing_table_values_rise[findMiddle(timing_table_transition_rise[0])]
            interp_rise = InterpolatedUnivariateSpline(xvals, fvals)
            cell_delay_rise.append(cell_delay_things(CellwithLoadbf[i].cell,interp_rise(CellwithLoadbf[i].load), 'none','none'))

        ###################For calculating Fall delay#######################
        for j in range(len(timing_table_capacitance_fall[0])):
            if (CellwithLoadbf[i].load == timing_table_capacitance_fall[0][j]):
                fall_delay = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])][j]

        if(fall_delay == -1):    #####if the capacitance was not found in the table, then we must do interpolation or extrapolation
            xvals_fall = timing_table_capacitance_fall[0]
            yvals_fall = [timing_table_transition_fall[0][2]]
            fvals_fall = timing_table_values_fall[findMiddle(timing_table_transition_fall[0])]
            interp_fall = InterpolatedUnivariateSpline(xvals_fall, fvals_fall)
            cell_delay_fall.append(cell_delay_things(CellwithLoadbf[i].cell,interp_fall(CellwithLoadbf[i].load), 'none','none'))
        rise_delay = -1
        fall_delay = -1


    for q in range(len(cell_delay_fall)):
        maxFallRise = max(cell_delay_fall[q].delay, cell_delay_rise[q].delay)
        cell_delay.append(cell_delay_things(cell_delay_fall[q].cell, maxFallRise, 'none','none'))

    for g in range(len(cell_delay)):
        if (cell_delay[g].delay > float(maxDelay)):
            cell_delay[g] = cell_delay[g]._replace(status = "Violating")
        else:
            cell_delay[g] = cell_delay[g]._replace(status = "Not Violating")

    #######################printing the recalculated delay##########################################
    with open('Bonus/Cell_delay_after_adding_buffers.txt', 'w') as f:
        for item7 in range(len(cell_delay)):
            f.write(str(cell_delay[item7].cell)+"     "+str(cell_delay[item7].delay)+  "\n")
####################################################################################################################