import linecache as lc
from numpy import random

lt=int(input("Enter the string length(even):"))

def conc(l,u,lis): #concatanates the sequence segments form different lines
    c=lis[0][l:50]
    for i in range(1,len(lis)-1):
        c=c+lis[i][:50]
    c=c+lis[len(lis)-1][:u]
    return c


# Python3 implementation to find the
# maximum number of occurrence of
# the overlapping count
 
# Function to find the maximum
# overlapping strings
def maxOverlap(S, T):
    max_overlap = 0
    i=0
    sm=0
    while(i< len(S)):
        j=0
        max_overlap=1
        while(j<len(T)):
            k=0
            while(i+k<len(S) and j+k< len(T) and T[j+k]==S[i+k]):
                k=k+1
            j=j+1
            if max_overlap<k:
                max_overlap=k
        if max_overlap>=5:
            sm=sm+max_overlap
        i=i+max_overlap
    return sm
    
        

def sequence(c,l,u):   #returns the parts of lines of interest as a list 
    l1=l//50 +2
    u1=u//50 +2
    lis=[]
    for i in range(l1,u1+1):
        line=lc.getline('chr/hg38'+c+'.fa',i)
        lis.append(line)
    line=""
    if len(lis)==1:
        line=lis[l%50:u%50+1]
    else:
        line=conc(l%50,u%50,lis)
    return line
    

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


f='ENCFF086CSF' #name of file without ".bed"
file1=open(f+'.bed','r')
#allowed is allowed starting indices of the sequences leaving leading and trailing Ns 
allowed=[(10100,248937600),(10100,242183500),(10100,198235500),(10100,19024500),(10100,181478250),(60100,170746000),(10100,159336000),(60100,145078500),(10100,138334750),(10100,133787400),(60100,135076600),(10100,133265300),(16000100,114354300),(18862750,106883700),(17000100,101981200),(10100,90228300),(489100,83247400),(10100,80263300),(60100,58607600),(60100,64334100),(60100,64334100),(5010050,46700000),(10510100,50808400),(10100,156030900),(10100,57217400)]
positive=[0]*24
ci='chr1'
i=0
si='ATGC'
lengthSum=[0]*24
numPos=[0]*24
count=0
sm=0
while count<10000:
    line=file1.readline()
    if not line:
        break
    #print(row)
    row=line.split('\t')
    c=row[0]
    li=int(row[1])
    ui=int(row[2])  #ui and li are starting and ending sequences as read from the bed file
    l=(ui+li)//2 - lt//2 
    u=(ui+li)//2 + lt//2  #l and u are the new starting and ending indices required length apart
    sc='ATGC'
    if c=='chrX':
        i=22
        s=u-l
        lengthSum[i]=lengthSum[i]+s
        numPos[i]=numPos[i]+1
        sc=sequence(c,l,u)
        m=maxOverlap(sc,si)
    elif c=='chrY':
        i=23
        s=u-l
        lengthSum[i]=lengthSum[i]+s
        numPos[i]=numPos[i]+1
        sc=sequence(c,l,u)
        m=maxOverlap(sc,si)
    elif RepresentsInt(c[3:]):
        i=int(c[3:])-1
        s=u-l
        lengthSum[i]=lengthSum[i]+s
        numPos[i]=numPos[i]+1
        sc=sequence(c,l,u)
        m=maxOverlap(sc,si)
    si=sc
    #print(m)
    sm=sm+m
    count=count+1
file1.close()

print(sm//100)
