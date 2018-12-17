import itertools
import igraph


#Combinatorial functions

def com_su(T,n):
    #Input: T subset of range(n).
    #Output: T's succesor under lexicographic ordering, range(k-1) (first element of [n]^(k-1)) if T is the maximum element.
    k = len(T)
    U = T[:]
    r = k
    while 0 < r and (T[r-1] == n-k+r):
        r -= 1
    if r == 0:
        return(range(1,k))
    for j in range(r,k+1):
        U[j-1] = (T[r-1] + j - r + 1)
    return(U)

def trans(T,S):
    #input: T a sublist of S.
    #output: the translation of T in [1,...,len(S)]^len(T)
    if S == range(1,len(S)+1):
        return(T)
    Y = []
    for t in T:
        Y.append(S.index(t)+1)
    return(Y)

def gen_comsu(T,V):
    #input: T subset(-list) of V.
    #output: T's succesor under lexicographic ordering for [V]^(k), k-1 if T is the maximum element
    S = [] #the succesor
    L = trans(T,V) #the set of indices of T in V
    n = len(V)
    for i in com_su(L,n):
        S.append(V[i-1])
    return(S)

def comblist(V,k):
    #Input: a set of vertices V and an integer k
    #the list of all the sublists of V of lenght k
    x = list(itertools.combinations(V,k))
    y = []
    for e in x:
        y.append(list(e))
    return(y)

def pr(S,R):
    #input: two list, S and R
    #output: true if S is contained in R, false otherwise
    for s in S:
        if not(s in R):
            return(False)
    return(True)

def co_suc(U,n):
    #Input: U subset of range(1,n+1) (it defines k = len(U)).
    #Output: T's succesor under colexicographic ordering, range(k-1) (first element of [n]^(k-1)) if T is the maximum element.
    T = sorted(U)
    T.reverse()
    k = len(U)
    i = 0
    while i < k-1 and T[k-1-i] == (T[k-i-2]-1):
        i += 1
    if i == k-1 and T[0] == n:
        return(range(1,k))
    T[k-1-i] += 1
    for j in range(i):
        T[k-1-j] = j + 1
    return(sorted(T))

def rank_co(T):
    #Input: T sublist of range(1,n+1)
    #Output: The rank of T under colexicographic ordering
    T = sorted(T)
    T.reverse()
    r = 0
    k = len(T)
    for i in range(k):
        r += binomial(T[i]-1,k-i)
    return(r)

def gen_co_suc(U,V):
    #Input: U sublist of V (define k = len(U))
    #Output: the succesor of U.
    S = [] #the succesor
    L = trans(U,V) #the set of indices of T in V
    n = len(V)
    for i in co_suc(L,n):
        S.append(V[i-1])
    return(S)

#starts regular zero forcing

def regular_colour_change(G,S):
    #Input:graph G, and a list S of black vertices.
    #Output: the set S2 obtained by applying one time the colour rule to G under S
    Y = []        #the list of white neighbors
    for v in S:        #iteratates over S to find a forcing vertex
        for w in G.neighbors(v):        #iterates over the neighbors of a black vertex to find it's white nighbors
            if not(w in S):
                Y.append(w)
                if len(Y) == 2:
                    Y = []
                    break
        if len(Y) == 1:
            S = S.append(Y[0])
            break
    return(S)

def posch(G,S):
    #Input: A graph G and a list S of black vertices.
    #Output: The set of vertices that are forced by S
    X = []
    can = []
    for s in S:
        for y in G.neighbors(s):
            if not(y in S) and not(y in can):
                X.append(y)
                if len(X) == 2:
                    break
        if len(X) == 1:
            can.append(X[0])
        X = []
    return(can)

def assoclab(G,S):
    #Input: A graph G and a set of black vertices S.
    #Output: The derived coloring of G
    x = list(S)[:]
    #print(S)#
    if len(posch(G,S)) != 0:
        x = x + posch(G,S)
        return(assoclab(G,x))
    else:
        return(sorted(S))

def firstzeroforc(G):
    #Input: A graph G
    #Output: First k lower or equal than n such that {v0,...,vk-1} is a zfs for G
    d = min(G.degree())
    for j in range(d,G.order()):
        if assoclab(G,G.vertices()[0:j]) == G.vertices():
            return(j)

def zeroforc2(G,a):
    #Input: A Graph G and an integer "a" such that $Z(G)\leq a$
    #Output: Z(G)
    x = [a]
    for l in range(a-1,0,-1):
        for s in itertools.combinations(G.vertices(),l):
            if assoclab(G,s) == G.vertices():
                x.append(l)
                break
        if not(l in x):
            return(l+1)

def zfs2(G,a):
    #Input: A graph G and an integer "a" greater than Z(G)
    #Output: ZFS for G of lenght a
    for s in itertools.combinations(G.vertices(),a):
        if assoclab(G,s) == G.vertices():
            return(sorted(list(s)))

def zero_forcing_num(G):
    #Input: A graph G.
    #Output: Z(G).
    print 'the zero forcing number of the graph is:'
    return(zeroforc2(G,firstzeroforc(G)))

def zfs(G):
    #Input: a Graph G
    #Output: a minimal ZFS set for G
    return(zfs2(G,zeroforc(G)))

def zeroforc1(G,a):
    #Input: A Graph G and an integer "a" such that $Z(G)\leq a.$
    #Output: Z(G).
    x = [a]
    V = G.vertices()
    zfs = V[:a]
    for m in range(a-1,0,-1):
        T = V[0:m]
        while len(T) == m:
            A = assoclab(G,T)
            if A == V:
                print('I did it for: ' + str(m))
                x.append(m)
                zfs = T
                break
            T = gen_co_suc(T,V)
        if not(m in x):
            print 'a zero forcing set is: '
            print zfs
            print 'zero forcing number is: '
            return(m+1)
    print 'the zero forcing number is:'
    return(1)

def zero_forcing(G):
    #Input: A graph G.
    #Output: Z(G).
    return(zeroforc1(G,firstzeroforc(G)))

def rg(t,n,p):
    #Input: $t\in {0,1}$, an integer n, $p\in [0,1]\cup\omega$
    #Output: generates a random graph of n vertices, depending on the value of t, if t = 0, p means the probability of an edge to be in the graph or if t = 1, p means the number of edges there's going to be in the graph#
    if t == 0:
        g = graphs.RandomGNP(n,p)
        return(g)
    if t == 1:
        g = graphs.RandomGNM(n,p)
        return(g)

#regular zero forcing stops here

def genall(n):
    #Input: An integer n
    #Ouput: Plots every graph (doesn't distinguish isomorphic graphs) on n vertices.
    g = Graph(n)
    z = comblist(n,2)
    for x in list((powerset(z))):
        g.add_edges(x)
        show(g)
        print
        g = Graph(n)

def listall(n):
    #Input: An integer n
    #Output: a list of all graphs (doesn't distinguish isomorphic graphs) on n vertices
    g = Graph(n)
    z = comblist(n,2)
    y = []
    for x in list((powerset(z))):
        g.add_edges(x)
        y.append(g)
        g = Graph(n)
    return(y)

#Z_+ starts here

def colrulplus(G,S):
    #Input: A graph G and a list of initial black vertices S.
    #Output: .
    H = G.copy()
    H.delete_vertices(S)
    can = []
    for V in H.connected_components():
        can = can + posch(G.subgraph(V+S),S)
    return(sorted(can))

def itcolplus(G,S):
    #Input: A graph G and a set of black vertices S.
    #Output: The associated labeling of S in G.
    s = S[:]
    if len(colrulplus(G,s)) != 0:
        s = s + colrulplus(G,s)
        return(itcolplus(G,s))
    else:
        return(sorted(s))

def firstplus(G):
    #Input: A graph G
    #Output: The first k such that {v1,...vk} forces plus V(G)
    d = min(G.degree(G.vertices()))
    for j in range(d,G.order()):
        if itcolplus(G,G.vertices()[0:j]) == G.vertices():
            return(j)

def zeroplus1(G,S):
    #Input: A Graph G and a positive zero forcing set S.
    #Output: Z_+(G).
    k = len(S)
    x = [k]
    V = G.vertices()
    zfs = S
    for l in range(k-1,0,-1):
        for s in itertools.combinations(V,l):
            s = list(s)
            if itcolplus(G,s) == V:
                x.append(l)
                zfs = s
                break
        if not(l in x):
            print('a minimal positive zero forcing set is: ' + str(zfs))
            return(l+1)
    return(1)

def zeroplus(G):
    #Input: A graph G
    #Output: Z_+(G)
    print 'The pos. zero forcing number of the graph is:'
    return(zeroplus1(G,G.vertices()[0:firstplus(G)]))

def Zplus_chordal(G):
    #Input: A Graph G.
    #Output:An optimal pos. zero forcing tree for G if G is chordal, an error otherwise.
    chord = G.is_chordal(True)
    if not(chord[0]):
        raise Exception("The graph is not chordal!")
    else:
        P = chord[1]
        B = []
        Bv = []
        R = []
        Wv = []
        V = G.vertices()
        n = len(V)
        if n == 1:
            return(set(set(V[0])))
        for i in range(n-1):
            Gi = G.subgraph(P[i:n])
            vi = P[i]
            Ni = Gi.neighbors(vi)
            Ei = Gi.edges_incident(vi,False)
            Ci = Gi.subgraph(Ni+[vi]).edges(False)
            for e in Ei:
                if not(e in R+B):
                    B.append(e)
                    for a in Ci:
                        if not(a in R+B):
                            R.append(a)
                    Wv.append(vi)
                    break
            if pr(Ei,R):
                Bv.append(vi)
        Bv.append(P[n-1])
        E = G.edges(False)
        for e in R:
            E.remove(e)
        T = G.subgraph(edges=E).connected_components()
        print('the tree covering and set of black vertices are: ')
        return(T,Bv)


G = graphs.CompleteGraph(5)
H = graphs.CompleteGraph(4)
F = G.cartesian_product(H)
F.plot()
#Zplus_chordal(F)
zeroforc(F)
zeroforc3(F)
zfs(F)
