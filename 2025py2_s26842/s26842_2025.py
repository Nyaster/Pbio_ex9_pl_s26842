import time,pylab as Y,pandas as D;from Bio import Entrez
p=print;e=input("E:");k=input("K:")or None;i=input("TID:")
n=input("MinL:");m=int(n)if n.isdigit()else None
x=input("MaxL:");M=int(x)if x.isdigit()else None
Entrez.email,Entrez.api_key,Entrez.tool=e,k,'M';W,Q,C,N=None,None,0,""
time.sleep(.34);h=Entrez.efetch(db="taxonomy",id=i,retmode="xml");R=Entrez.read(h);h.close();N=R[0]["ScientificName"];p(f"O:{N}({i})")
l=""
if m!=None and M!=None:l=f" AND {m}:{M if M>=m else'9'*10}[SLEN]"
elif m!=None:l=f" AND {m}:{'9'*10}[SLEN]"
elif M!=None:l=f" AND 1:{M}[SLEN]"
s=f"txid{i}[Organism]{l}";p(f"Q:{s}")
time.sleep(.34);h=Entrez.esearch(db="nucleotide",term=s,usehistory="y",idtype="acc");r=Entrez.read(h);h.close();C=int(r["Count"])
if C<1:p(f"No r {N}");exit()
p(f"{C}r");W,Q=r["WebEnv"],r["QueryKey"]
d=[];F=min(C,200);p(f"S:{F}")
time.sleep(.34);h=Entrez.esummary(db="nucleotide",webenv=W,query_key=Q,retstart=0,retmax=F);S=Entrez.read(h);h.close()
for e in S:d+=[{"a":e.get("AccessionVersion"),"l":int(e.get("Length",e.get("slen",0))),"d":e.get("Title")}]
if not d:p("No sum");exit()
b=i;D.DataFrame(d).to_csv(f"{b}r.csv",index=0);p(f"CSV:{b}r.csv");
f=D.DataFrame(d).sort_values("l", ascending=0)[:20]
Y.plot(f["a"],f["l"],"o-");Y.xlabel("A");Y.ylabel("L");Y.title(f"N:{len(f)}")
Y.xticks(rotation=90,fontsize=6);Y.tight_layout();Y.savefig(f"{b}p.png");Y.close();p(f"Plot:{b}p.png")
a=min(C,3);t=""
if a>0:p(f"GB:{a}");time.sleep(.34);h=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retmax=a,webenv=W,query_key=Q);t=h.read();h.close()
if t:
 with open(f"{b}s.gb","w")as z:z.write(t)