#!/usr/bin/env python
# coding: utf-8

# In[479]:


from liberty.parser import parse_liberty
liberty_file = "/Users/ayashaker/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())
print(library)


# In[480]:


fname = "/Users/ayashaker/Desktop/rca4.rtlnopwr.v"


# In[481]:


#####Parsing Netlist##########
import re


# In[482]:


moduledefintion = ""
inn = []
out = []
wire = []
str = ""
functions = []
for code in open(fname,"r"):
    line = code
    if("input" in line):
        line = re.sub("input", '', line) 
        inn.append(line)
    elif("out" in line):
        line = re.sub("output", '', line) 
        out.append(line)
    elif("wire" in line):
        line = re.sub("wire", '', line) 
        wire.append(line)
    elif("endmodule" in line):
        str = "end of module"
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
    del n[1]
    dictList.append(n)
    


# In[484]:


for i in range(len(dictList)):
    for j in range(2):
        n =  re.sub("\n", '', dictList[i][j])
        dictList[i][j] = n
        n =  re.sub(";", '', dictList[i][j])
        dictList[i][j] = n


# In[485]:


for i in range(len(dictList)):
    n = dictList[i][1].split()
    dictList[i][1] = n
for i in range(len(dictList)):
    del dictList[i][1][0]
    del dictList[i][1][len(dictList[i][1])-1]
        


# In[486]:


dictList


# In[ ]:




