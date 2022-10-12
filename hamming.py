from numpy import random
import linecache as lc
def hamming(s1,s2):
    h=0
    s1=s1.upper()
    s2=s2.upper()
    for i in range (0,len(s1)):
        if s1[i]!=s2[i]:
            h=h+1
    return h

def conc(l,u,lis):
    c=lis[0][l:50]
    for i in range(1,len(lis)-2):
        c=c+lis[i][:50]
    c=c+lis[len(lis)-1][:u]
    return c

s1='ggttccaagttagactttgccactgagagaaactcactcaaggcttggaatgtgggagtgaagcagaagtcattgttctcaggaggctgtgacagtcagaTAGTGTGGCttcacagatctctctatgacgtccacttccgattagccttatgcaactagacctggtgacttctcagaattcctgtgagcttttgtttcctcacctcctccagtgtttcaagttaaagtcattagagacatttgtgtgatcctccagcatgaatcttccagactttcactctctagctcccactttttatg'


x=random.randint(48937071,155084550,1000)
lines=[]
for l in x:
    u=l+200
    l1=l//50 +2
    u1=u//50 +2
    lis=[]
    for i in range(l1,u1+1):
        lis.append(lc.getline('chr/hg38chrX.fa',i))
    line=""
    if len(lis)==1:
        line=lis[l%50:u%50+1]
    else:
        line=conc(l%50,u%50,lis)
    lines.append(line)
for i in range(0,len(lines)-1):
    print(hamming(lines[i],s1))

