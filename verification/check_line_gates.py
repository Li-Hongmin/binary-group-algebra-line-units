#!/usr/bin/env python3
"""Independent brute-force checks for 'Line Gates and Orbit Bounds in Binary
Group Algebras'. Standard library only. Numerical/finite checks, not proofs."""
import itertools, math, sys
from fractions import Fraction as Fr
POLY={2:0b111,3:0b1011,4:0b10011}
def gf(n):
    q=1<<n; P=POLY[n]
    def mul(a,b):
        r=0
        while b:
            if b&1: r^=a
            b>>=1; a<<=1
            if a&q: a^=P
        return r
    M=[[mul(a,b) for b in range(q)] for a in range(q)]
    inv=[0]*q
    for a in range(1,q):
        for b in range(1,q):
            if M[a][b]==1: inv[a]=b
    return q,M,inv
FAIL=[]
def check(cond,msg):
    if not cond: FAIL.append(msg); print("FAIL",msg)
# points of F_q^rho encoded as ints, n bits per coordinate; addition = XOR
def enc(v,n): return sum(c<<(n*i) for i,c in enumerate(v))
def dec(x,n,rho): return tuple((x>>(n*i))&((1<<n)-1) for i in range(rho))
def smul(g,x,n,rho,M): return enc([M[g][c] for c in dec(x,n,rho)],n)
def line(v,n,rho,M): return frozenset(smul(g,enc(v,n),n,rho,M) for g in range(1<<n))
def prod(A,B):
    r=set()
    for a in A:
        for b in B: r^={a^b}
    return frozenset(r)
def units(Lpts):
    L=sorted(Lpts)
    for bits in range(1<<len(L)):
        if bin(bits).count('1')%2: yield frozenset(L[i] for i in range(len(L)) if bits>>i&1)

# ---- Sections 3-4: gate, minimal support, survival laws (rho=2, W=G) ----
def gate_checks(n):
    q,M,inv=gf(n); rho=2
    lines=[line((1,m),n,rho,M) for m in range(q)]+[line((0,1),n,rho,M)]
    one=frozenset([0])
    for (i1,i2,i3) in itertools.combinations(range(len(lines)),3):
        L1,L2,L3=lines[i1],lines[i2],lines[i3]
        # find a in L1, b in L2 with a+b in L3 (a,b nonzero)
        a,b=next((a,b) for a in L1 if a for b in L2 if b and (a^b) in L3)
        H=frozenset([0,a,b,a^b]); K=H
        D=one
        for L in (L1,L2,L3): D=prod(D,one^L)
        check(prod(K,D)==K,f"KD=K q={q}")
        if (i1,i2,i3)!=(0,1,2) and q>4: continue
        # minimal support: no nonzero B with support<=3 fixed by D
        pts=range(q*q)
        for k in (1,2,3):
            for S in itertools.combinations(pts,k):
                if 0 not in S: continue  # translation invariance
                B=frozenset(S)
                check(prod(B,D)!=B,f"support {k} fixed q={q}")
        # survival laws: one and two units on every line pair, three on defining lines
        h=lambda L: len(H&L)
        U={id(L):list(units(L)) for L in lines}
        for L in lines:
            cnt=sum(prod(K,A)==K for A in U[id(L)])
            check(Fr(cnt,len(U[id(L)]))==Fr(2)**(1-q//h(L)),f"one-unit q={q}")
        for La,Lb in itertools.combinations(lines,2):
            KA={}
            for A in U[id(La)]:
                x=prod(K,A); KA[x]=KA.get(x,0)+1
            cnt=sum(c for x,c in KA.items() for B in U[id(Lb)] if prod(x,B)==K)
            pred=Fr(2)**(1+Fr(4,h(La)*h(Lb))-Fr(q,h(La))-Fr(q,h(Lb)))
            check(Fr(cnt,len(U[id(La)])*len(U[id(Lb)]))==pred,f"two-unit q={q}")
        d=q//2
        def law(seq):
            dist={K:1}
            for L in seq:
                nd={}
                for x,c in dist.items():
                    for A in U[id(L)]:
                        y=prod(x,A); nd[y]=nd.get(y,0)+c
                dist=nd
            tot=1
            for L in seq: tot*=len(U[id(L)])
            return Fr(dist.get(K,0),tot)
        seqs=[(L1,),(L1,L1),(L1,L2),(L1,L1,L2),(L1,L2,L3)]
        if q==4: seqs+=[(L1,L2,L3,L1),(L1,L2,L2,L3,L3)]
        for seq in seqs:
            j=len(set(map(id,seq)))
            pred={1:Fr(2)**(1-d),2:Fr(2)**(2-2*d),3:d*Fr(2)**(-3*(d-1))}[j]
            check(law(seq)==pred,f"three-line law q={q} t={len(seq)} j={j}")
        # Proposition 4.4: arbitrary line families; |P| = number of reachable values K*prod
        def reach(seq):
            R={K}
            for L in seq: R={prod(x,A) for x in R for A in U[id(L)]}
            return R
        fams=[s3 for s3 in itertools.combinations(lines,3)] if q==4 else [(lines[0],lines[3],lines[5]),(lines[1],lines[2],lines[6])]
        fams=list(fams)+[(lines[0],lines[1],lines[2],lines[3])]+[(lines[0],lines[3],lines[0])]
        for seq in fams:
            P=len(reach(seq)); pr=law(seq); single=Fr(1)
            for L in seq: single*=Fr(2)**(1-q//h(L))
            check(pr==Fr(1,P),f"family law 1/|P| q={q}")
            check(P&(P-1)==0,f"family |P| power of 2 q={q}")
            check(pr>=single,f"family >= product q={q}")
            check(law(tuple(dict.fromkeys(seq)))==pr,f"family repetition q={q}")
    print(f"gate/survival checks done q={q}")

# ---- Section 5: exposed columns, exact formula and orbit bound ----
def rank_rows(rows):
    r=0; rows=list(rows); basis=[]
    for v in rows:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    return len(basis)
def exposure_checks(n,rho,maxE=None):
    q,M,inv=gf(n)
    worst=0
    for dvec in itertools.product(range(1,q),repeat=rho-1):
        if all(x==1 for x in dvec): continue
        J=[j for j in range(rho-1) if dvec[j]!=1]; s=len(J); t=rho-1-s; p=Fr(1,2**(rho-1))
        ell=lambda j,a,x: M[1^dvec[j]][a]^M[dvec[j]][x]
        nX=(rho-1)*q
        # H_X(a,x) via the direct definition with x=a+v_rho g, v_rho=1
        Hcache=[]
        for X in range(1<<nX):
            bit=lambda j,u: (X>>(j*q+u))&1
            H=[[0 if x==a else int(all(bit(j,ell(j,a,x)) for j in range(rho-1))) for x in range(q)] for a in range(q)]
            Hcache.append(H)
        Es=[E for Msz in range(1,q//2+1) for E in itertools.combinations(range(q),Msz)]
        if maxE: Es=Es[:maxE]
        for E in Es:
            Msz=len(E)
            # greedy exposed set (increasing label)
            nbr={x:set() for x in range(q)}
            for j in J:
                c0=M[inv[dvec[j]]][1^dvec[j]]
                for a,b in itertools.permutations(E,2):
                    for x in range(q):
                        y=x^M[c0][a^b]
                        if y!=x: nbr[x].add(y); nbr[y].add(x)
            F=[]; banned=set(E)
            for x in range(q):
                if x not in banned: F.append(x); banned|=nbr[x]|{x}
            RM=(q-Msz)//(s*Msz*(Msz-1)+1)
            check(len(F)>=RM,"R_M")
            for cbits in range(1,1<<Msz):
                C=[E[i] for i in range(Msz) if cbits>>i&1]; m=len(C)
                cnt=sum(all(sum(H[a][x] for a in C)%2==0 for x in F) for H in Hcache)
                prob=Fr(cnt,len(Hcache))
                pred=(1-Fr(1,2**t)*(1-(1-Fr(2,2**s))**m)/2)**len(F)
                check(prob==pred,f"exact exposure q={q} rho={rho} d={dvec} E={E} c={C}")
                check(prob<=(1-p)**len(F),"exposure bound")
            # rank-deficiency and orbit bound
            rd=sum(rank_rows(sum(H[a][x]<<x for x in range(q)) for a in E)<Msz for H in Hcache)
            check(Fr(rd,len(Hcache))<=(2**Msz-1)*(1-p)**RM,"rank bound")
        # exact orbit-constancy probability for the full target line
        tot=Fr(0)
        for H in Hcache:
            rows=[sum(H[a][x]<<x for x in range(q)) for a in range(q)]
            r=rank_rows(rows)
            # #Y with HY=0 is 2^{q-r}; HY=all-ones solvable iff ones in column space
            cols=[sum(H[a][x]<<a for a in range(q)) for x in range(q)]
            r2=rank_rows(cols+[(1<<q)-1])
            tot+=Fr(2**(q-r)*(1 if r2>r else 2),2**q)
        P=tot/len(Hcache)
        best=min(2**Msz*math.exp(-float(p)*((q-Msz)//(s*Msz*(Msz-1)+1)))+2**(1-Msz) for Msz in range(1,q//2+1))
        check(float(P)<=best+1e-12,f"orbit upper q={q}")
        worst=max(worst,float(P)/best)
    print(f"exposure checks done q={q} rho={rho} (max ratio orbit/bound {worst:.3f})")

# ---- Section 6: repeated-slope family ----
def poisson_checks(n):
    q,M,inv=gf(n)
    for gam in range(2,q):
        # rank facts over all active sets A
        for Abits in range(1<<q):
            A=[h for h in range(q) if Abits>>h&1]
            Hrows=[sum(1<<x for h in A if h!=a for x in [a^M[gam][h^a]]) for a in range(q)]
            # direct definition check: z(a)=sum_{h in A,h!=a} Y_{a+gamma(h-a)}
            Brows=[sum(1<<x for x in range(q) if (M[inv[gam]][x]^M[1^inv[gam]][a]) in A) for a in range(q)]
            check(all((Hrows[a]^Brows[a])==((1<<a) if a in A else 0) for a in range(q)),"diag")
            eps=len(A)%2; l=0
            for h in A: l^=h
            rt=rank_rows([sum(1<<(g^h) for h in A) for g in range(q)])
            if eps: check(rt==q,"odd rank")
            elif l: check(rt>=q//2,"rank q/2")
            check(rank_rows(Hrows)>=rt-len(A),"rank ineq")
        if n>3: continue
        for rho in (3,4,5):
            p=Fr(1,2**(rho-1)); P=Fr(0); J=Fr(0)
            for Abits in range(1<<q):
                A=[h for h in range(q) if Abits>>h&1]
                w=p**len(A)*(1-p)**(q-len(A))
                rows=[sum(1<<x for h in A if h!=a for x in [a^M[gam][h^a]]) for a in range(q)]
                r=rank_rows(rows)
                cols=[sum((rows[a]>>x&1)<<a for a in range(q)) for x in range(q)]
                r2=rank_rows(cols+[(1<<q)-1])
                P+=w*Fr(2**(q-r)*(1 if r2>r else 2),2**q)
                l=0
                for h in A: l^=h
                if len(A)%2==0 and l==0: J+=w
            check(J==(1+(1-2*p)**q+(2*q-2)*(1-2*p)**(q//2))/(2*q),"jet formula")
            bb=(1-2*p)**(q//2)
            check(J==bb+(1-bb)**2/(2*q),"jet rewrite")
            check(J<=(1-p)**q+Fr(1,2*q),"jet uniform bound")
            check((1-p)**q<=P,"poisson lower")
            ub=float((1-p)**q)+1/(2*q)+math.exp(-q/72)+2**(1-q/6)
            check(float(P)<=ub,"poisson upper")
    print(f"repeated-slope checks done q={q}")

# ---- Section 7: association lower bounds (q=4, rho=2,3) ----
def association_checks(n,rho):
    q,M,inv=gf(n); theta=Fr(1,2**rho)
    # gate inside first two coordinates
    base=[line((1,m)+(0,)*(rho-2),n,rho,M) for m in range(q)]+[line((0,1)+(0,)*(rho-2),n,rho,M)]
    L1,L2,L3=base[:3]
    a,b=next((a,b) for a in L1 if a for b in L2 if b and (a^b) in L3)
    H=frozenset([0,a,b,a^b]); K=H
    dirs=set()
    for v in itertools.product(range(q),repeat=rho):
        if any(v): dirs.add(line(v,n,rho,M))
    dirs=list(dirs)
    nX=rho*q
    def Qval(B,X):
        return sum(all((X>>(j*q+dec(g,n,rho)[j]))&1 for j in range(rho)) for g in B)%2
    def const_prob(B,T):
        orbit=[frozenset(t^g for g in B) for t in T]
        return Fr(sum(len({Qval(Bt,X) for Bt in orbit})==1 for X in range(1<<nX)),1<<nX)
    def const_event(B,T,X):
        return len({Qval(frozenset(t^g for g in B),X) for t in T})==1
    Ms=dirs if rho==2 else dirs[:6]
    for T in (dirs if rho==2 else dirs[:4]):
        for Mline in Ms:
            avg=Fr(0); U=list(units(Mline))
            for A in U:
                B=prod(K,A); P=const_prob(B,T); avg+=P/len(U)
                check(P>=(1-theta)**(4*q*len(A)),"fixed lower")
                if Mline==T:
                    check(P==const_prob(K,T),"same-line equality")
                    check(all(const_event(B,T,X)==const_event(K,T,X) for X in range(1<<nX)),"same-line equality of events")
                    check(P>=(1-theta)**(4*q//len(H&T)),"same-line lower")
            r=(1-theta)**(4*q)
            check(avg>=((1+r)**q-(1-r)**q)/2**q,"average lower")
    print(f"association checks done q={q} rho={rho}")

if __name__=="__main__":
    gate_checks(2); gate_checks(3)
    exposure_checks(2,2); exposure_checks(2,3); exposure_checks(3,2)
    poisson_checks(3); poisson_checks(4)
    association_checks(2,2); association_checks(2,3)
    print("ALL CHECKS PASSED" if not FAIL else f"{len(FAIL)} FAILURES"); sys.exit(1 if FAIL else 0)
