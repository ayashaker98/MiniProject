#!/usr/bin/env python
# coding: utf-8

# In[15]:


########Parsing Lib File#######
from liberty.parser import parse_liberty
liberty_file = "/Users/ayashaker/Desktop/osu035.lib"
library = parse_liberty(open(liberty_file).read())
print(str(library))


# In[312]:


fname = "/Users/ayashaker/Desktop/rca4.rtlnopwr.v"


# In[346]:


#####Parsing Netlist##########
import re


# In[347]:


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


# In[350]:


dictList = []
for code in functions:
    n = code.split(maxsplit=2)
    del n[1]
    dictList.append(n)
    


# In[351]:


for i in range(len(dictList)):
    for j in range(2):
        n =  re.sub("\n", '', dictList[i][j])
        dictList[i][j] = n
        n =  re.sub(";", '', dictList[i][j])
        dictList[i][j] = n


# In[352]:


for i in range(len(dictList)):
    n =  re.sub("\n", '', dictList[i][1])
    dictList[i][1] = n

        


# In[ ]:


dictList


# In[ ]:




