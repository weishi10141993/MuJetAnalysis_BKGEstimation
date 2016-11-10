import math
from math import *
import subprocess

#Fuction that retirn true is a word is contained in file1 ina  given column
def check(fileMine,wordsMine,column):
  with open(fileMine) as f1:
    for line in f1:
      words=line.split()
      if(str(words[column])==str(wordsMine)):
        return True
    return False

#Files containing for each histograms a line "row run event lumi"
file1="scan1" 
file2="scan2"
column=2 #column on the event number
numDup=0
veto="e+09"
InitialLine=0

#Open file2 that in genral is shorter
with open(file1) as f:
  for line in f:
    words=line.split()
    if( veto not in words[column] ):
      InitialLine=InitialLine+1
      matched = check(file2,str(words[column]),column)
      if(matched):
        print str(words[column]) + " is duplicate!"
        numDup=numDup+1
print "You have " + str(numDup) + " duplicate on " + str(InitialLine) + " lines"
