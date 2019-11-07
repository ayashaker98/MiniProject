#!/usr/bin/env python
# coding: utf-8

# In[19]:


from liberty.parser import parse_liberty
from liberty.types import Group,select_cell,select_pin,select_timing_table
liberty_file = "/Users/ayashaker/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())


# In[20]:


fname = "/Users/ayashaker/Desktop/rca4.rtlnopwr.v"


# In[21]:


#####Parsing Netlist##########
import re


# In[22]:


maxFanout = input("Please Enter the maximum Fanout for any cell : ") 


# In[23]:


moduledefintion = ""
inn = []
out = []
wire = []
strr = ""
functions = []
for code in open(fname,"r"):
    line = code
    if("input" in line):
        line = re.sub("input", '', line) 
        inn.append(line)
        line = ""
    elif("output" in line):
        line = re.sub("output", '', line) 
        out.append(line)
        line = ""
    elif("wire" in line):
        line = re.sub("wire", '', line) 
        wire.append(line)
        line = ""
    elif("endmodule" in line):
        strr = "end of module"
        line = ""
    elif ("module" in line):
        moduledefintion = line
        line = ""
    elif(line == "\n"):
        donothing = ""
        line = ""
    else:
        functions.append(line)
        line = ""

NumberofCells = len(functions)


# In[24]:


functions


# In[25]:


dictList = []
for code in functions:
    n = code.split(maxsplit=2)
    del n[0]
    dictList.append(n)
    


# In[26]:


for i in range(len(dictList)):
    for j in range(2):
        n =  re.sub("\n", '', dictList[i][j])
        dictList[i][j] = n
        n =  re.sub(";", '', dictList[i][j])
        dictList[i][j] = n


# In[27]:


for i in range(len(dictList)):
    n = dictList[i][1].split()
    dictList[i][1] = n
for i in range(len(dictList)):
    del dictList[i][1][0]
    del dictList[i][1][len(dictList[i][1])-1]


# In[28]:


dict2 = dict(dictList)


# In[ ]:





# In[29]:


from collections import namedtuple


# In[33]:


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

Fanoutnumber = 0
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
        Fanoutnumber = 0

for g in range(len(Fanout)):
    if (Fanout[g].fanout > int(max(Fanout))):
        Fanout[g] = Fanout[g]._replace(status = "Violating")
    else:
        Fanout[g] = Fanout[g]._replace(status = "Not Violating")

# print(Fanout)
# print(len(Fanout))    
# print(allpins[6].pin.split('_')[0])
cellgroup = library.get_group('cell', allpins[5].cell.split('_')[0])
pin = cellgroup.get_group('pin', allpins[5].pin)
# print(allpins[5].pin)
timing_table = select_timing_table(pin, allpins[4].pin, 'cell_fall')
# print(timing_table)
timing_table_transition = timing_table.get_array("index_1")
timing_table_capacitance = timing_table.get_array("index_2")
timing_table_values = timing_table.get_array("values")

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
    CapPins.clear()

# print(timing_table_transition)
# print(timing_table_capacitance)
# print(timing_table_values)


# In[100]:


#!/usr/bin/env python
# coding: utf-8

# In[479]:

from liberty.types import Group, select_cell, select_pin, select_timing_table
from liberty.parser import parse_liberty
from collections import namedtuple
import array
from scipy.interpolate import interp1d
liberty_file = "/Users/ayashaker/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())
# print(library)

fname = "/Users/ayashaker/Desktop/rca4.rtlnopwr.v"

import re

maxFanout = 4 #input("Please Enter the maximum Fanout for any cell : ") 

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

Fanoutnumber = 0
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
        Fanoutnumber = 0

for g in range(len(Fanout)):
    if (Fanout[g].fanout > int(maxFanout)):
        Fanout[g] = Fanout[g]._replace(status = "Violating")
    else:
        Fanout[g] = Fanout[g]._replace(status = "Not Violating")

# print(Fanout)
# print(len(Fanout))    
# print(allpins[6].pin.split('_')[0])
cellgroup = library.get_group('cell', allpins[5].cell.split('_')[0])
pin = cellgroup.get_group('pin', allpins[5].pin)
# print(allpins[5].pin)
timing_table = select_timing_table(pin, allpins[4].pin, 'cell_fall')
# print(timing_table)
timing_table_transition = timing_table.get_array("index_1")
timing_table_capacitance = timing_table.get_array("index_2")
timing_table_values = timing_table.get_array("values")

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

# print(timing_table_transition)
# print(timing_table_capacitance)
# print(timing_table_values)


# In[45]:


print(Fanout)


# In[43]:


print(functions)


# In[46]:


FanoutCouple 


# In[137]:


import collections
Fanout = collections.namedtuple('Fanoutpairs', 'gate1 gate2 outputpin1 inputpin2 pin1 pin2 inputcap')
Fanoutpairs = []
for i in range(len(allpins)):
    for j in range(len(allpins)):
        if ((allpins[i].parameter  == allpins[j].parameter)&(allpins[i].direction == 'output')&(allpins[j].direction =='input')):
            print("Parameter: " , allpins[i].parameter)
            print("Gate 1 : ",allpins[i].cell)
            print("Gate 2 : ",allpins[j].cell)
            Fanoutpairs.append(Fanout(allpins[i].cell,allpins[j].cell,allpins[i].parameter,allpins[j].parameter,allpins[i].pin,allpins[j].pin,0))


# In[113]:


CapPins[0].cellname == Fanoutpairs[0].gate2 && CapPins[0].pinwithmaxcap == Fanoutpairs[0].pin2
    CapPins[0].inputcap = CapPins[0].inputcap._replace(Fanoutpairs[0].capacitance)


# In[144]:


for i in range(len(Fanoutpairs)):
    for j in range(len(CapPins)):
        if((CapPins[j].cellname == Fanoutpairs[i].gate2) & (CapPins[j].pinwithmaxcap == Fanoutpairs[i].pin2)):
            Fanoutpairs[i]= Fanoutpairs[i]._replace( inputcap = CapPins[j].capacitance)


# In[147]:


Fanoutpairs


# In[149]:


CapPins


# In[ ]:




