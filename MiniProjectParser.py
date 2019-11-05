#!/usr/bin/env python
# coding: utf-8

# In[479]:

from liberty.types import Group, select_cell, select_pin, select_timing_table
from liberty.parser import parse_liberty
from collections import namedtuple
import array
liberty_file = "/Users/noha/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())
# print(library)

fname = "/Users/noha/Desktop/rca4.rtlnopwr.v"

import re

maxFanout = input("Please Enter the maximum Fanout for any cell : ") 

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
Fanout = []
for l in range(len(allpins)):
    if (allpins[l].direction=="output"):
        count = count + 1
        for h in range(len(allpins)):
            if (allpins[h].parameter == allpins[l].parameter):
                Fanoutnumber = Fanoutnumber + 1
        Fanout.append(Fanoutnumber)
        Fanoutnumber = 0

print(Fanout)
print(len(Fanout))

# for g in range(len(Fanout)):
#     if (Fanout[g] > maxFanout):


# print(select_cell(library, allpins[0].cell))
