def hilbert_function(I):
    '''
    Given an ideal, compute the hilbert function
    '''
    B=I.normal_basis()              #B is the basis of R/I as a vector space
    D=[]                            #D will record the degree of each element of B
    for i in range(len(B)):
        D.append(B[i].degree())
    H=[]
    for j in range(len(B)):        #H will record the vector space dimension in
        H.append(D.count(j))       #each degree. At the end, zeros are dropped.
    while 0 in H:
        H.remove(0)
    return H

def Hilbert_Function(I):
    B=I.normal_basis()              #B is the basis of R/I as a vector space
    D=[]                            #D will record the degree of each element of B
    for i in range(len(B)):
        D.append(B[i].degree())
    H=[]
    for j in range(len(B)):        #H will record the vector space dimension in
        H.append(D.count(j))       #each degree. At the end, zeros are dropped.
    while 0 in H:
        H.remove(0)
    return H



def matrix_rep(I,f):
    '''
    returns the matrix representation of a given function f with ideal I
    '''
    basis = I.normal_basis()
    length = len(basis)
    L = []
    for term in basis:
        l = []
        mod_term = f*term.reduce(I)
        if mod_term != 0:
            mod_term_list = mod_term.monomials()
            mod_term_co = mod_term.coefficients()
            for term in basis:
                if term in mod_term_list:
                    co = mod_term_co[mod_term_list.index(term)]
                    l.append(co*1)
                else:
                    l.append(0)
        else:
            l = [0]*length
        L.append(l)
    M = matrix(L)
    return M.transpose()



def jordan_form(I,f):
    '''
    returns the jordan form of the matrix representation of a given function f and ideal I
    '''
    M = matrix_rep(f,I)
    return M.jordan_form()


def diff(L):
    '''
    Return a list of the differences between the elements in a list L
    '''
    D=[]
    for i in range(len(L)-1):
        D.append(L[i+1]-L[i])
    D.append(0)
    return(D)


def Jordan_Type(I,l):
    '''
    Return the jordan type given an ideal I and an l
    '''
    n=len(I.normal_basis())
    corank=[0]
    i=1
    while corank[-1]<n:
        f=l^i
        J=I+ideal(f)
        corank.append(len(J.normal_basis()))
        i=i+1
    D=diff(corank)
    D=Partition(D)
    Q=D.conjugate()
    return(Q)


def Jordan_Seq(I,l):
    n=len(I.normal_basis())
    rankseq=[]
    i=1
    while len(rankseq)<n:
        f=l^i
        J=I+ideal(f)
        rankseq.append(n-len(J.normal_basis()))
        i=i+1
    return rankseq

def Jordan_Type2(I,l):
    L=Jordan_Seq(I,l)
    j=len(I.normal_basis())-L[0]
    Jd=[j]
    i=1
    while i < len(L):
        j=L[i-1]-L[i]
        Jd.append(j)
        i=i+1
    J=Partition(Jd).conjugate()
    return J



'''
Hilbert function data collection
'''


def hilbert_function_I(m,n):
    '''
    Returns a list of hilbert functions for ideals I = (x^i,y^j) where 0 < i < m and 0 < j < n,
    '''
    R.<x,y> = PolynomialRing(QQbar,2)
    L = []
    for j in range(2,n):
        for i in range(2,m):
            I = ideal(x^i,y^j)
            L.append(hilbert_function(I))
    return L

def hilbert_function_J(n):
    '''
    Returns a list of hilbert functions for ideals J = (x^i,y^(i+1),xy) where 0 < i < n
    '''
    R.<x,y> = PolynomialRing(QQbar,2)
    L = []
    for i in range(2,n):
        I = ideal(x^i,y^(i+1),x*y)
        L.append(hilbert_function(I))
    return L

def hilbert_function_K(n):
    '''
    Returns a list of hilbert functions for ideals K = (x^i - y^i, x*y) where 0 < i < n
    '''
    R.<x,y> = PolynomialRing(QQbar,2)
    L = []
    for i in range(2,n):
        I = ideal(x^i - y^i,x*y)
        L.append(hilbert_function(I))
    return L


def hilbert_function_L(n):
    '''
    Returns a list of hilbert functions for ideals L = (x^2 - xy, y^i) where 0 < i < n
    '''
    R.<x,y> = PolynomialRing(QQbar,2)
    L = []
    for i in range(2,n):
        I = ideal(x^2 - x*y, y^i)
        L.append(hilbert_function(I))
    return L

def hilbert_function_L2(n):
    '''
    Returns a list of hilbert functions for ideals L = (x^3 - xy, y^i) where 0 < i < n
    '''
    R.<x,y> = PolynomialRing(QQbar,2)
    L = []
    for i in range(2,n):
        I = ideal(x^3 - x*y, y^i)
        L.append(hilbert_function(I))
    return L

def hilbert_function_L3(n,m):
    '''
    Returns a list of hilbert functions for ideals L = (x^3*y^3, x^i, y^j) where 0 < i < n
    '''
    R.<x,y> = PolynomialRing(QQbar,2)
    L = []
    for j in range(2,n):
        for i in range(2,m):
            I = ideal(x^3*y^3, x^i, y^j)
            L.append(hilbert_function(I))
    return L

def let_b_to_0(I,k):
    '''
    Returns a list of jordan types for a singular ideal with l = x + by as b --> 0 (approximated by a limit 1/k)

    conclusion -- NOTHING CHANGES!! The coefficient changes doesn't influence the jordan type
    '''
    b = 1
    L = []
    while b > 1/k:
        L.append(Jordan_Type(I, x + b*y))
        b = b*(1/2)
    return L

def hilbert_function_Ji_2(n):
    '''
    Returns a list of hilbert functions for ideals J = (x^i, y^i, z^i, xy, xz, zy) where 2 < i < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    for i in range(2,n):
        I = ideal(x^i, y^i, z^i, x*y, x*z, z*y)
        print(Hilbert_Function(I))
        L.append(Hilbert_Function(I))
    return L

def hilbert_function_Ji_3(n):
    '''
    Returns a list of hilbert functions for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Hilbert_Function(I))
        L.append(Hilbert_Function(I))
    return L

def hilbert_function_N(n):
    '''
    Returns a list of hilbert functions for ideals
     N = (x^i, y^j - x^(j-1)*z, z^2, x^(i-j+1)*y^(j-k), y^(j-k)*z) where k < j <= i
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    for k in range(n):
        for j in range(k + 1,n):
            for i in range(j,n):
                I = ideal(x^i, y^j - x^(j-1)*z, z^2, x^(i-j+1)*y^(j-k), y^(j-k)*z)
                print(i)
                print(Hilbert_Function(I))
                L.append(Hilbert_Function(I))
    return L


def jordan_types_Ji_3_x(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = x
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    l = x
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None

def jordan_types_Ji_3_y(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = y
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    l = y
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None

def jordan_types_Ji_3_z(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = z
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    l = z
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None

def jordan_types_Ji_3_xy(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = x + y
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    l = x + y
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None

def jordan_types_Ji_3_xz(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = x + z
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    l = x + z
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None

def jordan_types_Ji_3_yz(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = y + z
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    L = []
    l = y + z
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None

def jordan_types_Ji_3_xyz(n):
    '''
    Returns a list of jordan types for ideals J = (x^i, y^i, z^i, xyz) where 3 < i < n
    for l = x + y + z
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    l = x + y + z
    print('l=', l)
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        print(Jordan_Type2(I,l))

    return None


def possible_jordan_types(H):
    H.sort()
    H.reverse()
    Dk = Partition(H).conjugate()
    N = sum(H)
    L = []
    for i in Partitions(N):
        if Dk.dominates(i) and sum(i) != len(i):
            L.append(i)
    return L



def perfect_pairing(I,f,L):
    P=[]
    for i in range(len(L)):
        if f*L[i] not in I:
            P.append(i)
    return(P)


def is_AG(I):
        if I.dimension()!=0:
            return False
        else:
            H=Hilbert_Function(I)
            socle_degree=len(H)-1
            n=floor(socle_degree/2)
            for i in range(n+1):
                B=I.normal_basis(i)
                C=I.normal_basis(socle_degree-i)
                if len(B) != len(C):
                    return False #non-symmetric Hilb
                else:
                    for i in range(n+1):
                        B=I.normal_basis(i)
                        C=I.normal_basis(socle_degree-i)
                        for j in range(len(B)):
                            if len(perfect_pairing(I,B[j],C)) == 0:
                                #print('No perfect pairing in degree', i)
                                return False #not perfect pairing

                    return True

def is_AG2(I):
        if I.dimension()!=0:
            return('The algebra is not Artinian.')
        else:
            H=Hilbert_Function(I)
            socle_degree=len(H)-1
            n=floor(socle_degree/2)
            for i in range(n+1):
                B=I.normal_basis(i)
                C=I.normal_basis(socle_degree-i)
                if len(B) != len(C):
                    return('The algebra is not Gorenstein.') #non-symmetric Hilb
                else:
                    for i in range(n+1):
                        B=I.normal_basis(i)
                        C=I.normal_basis(socle_degree-i)
                        for j in range(len(B)):
                            if len(perfect_pairing(I,B[j],C)) == 0:
                                print('No perfect pairing in degree', i)
                                return('The algebra is not Gorenstein.')#not perfect pairing

                    return('The algebra is Gorenstein.')


def J_i_2_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    J = (x^i, y^i, z^i, xy, xz, zy) for 3 <= i < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,y+z,x+z,x+y+z]
    for i in range(3,n):
        I = ideal(x^i, y^i, z^i, x*y, x*z, z*y)
        L = get_list_info(I,P)
        I1 = I.gens()
        D[I1] = L
        print((I1),L[0])
        for list in L[1]:
            print(list)
        print(L[2])
    return D

def J_i_3_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    J = (x^i, y^i, z^i, xyz) for 3 <= i < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,y+z,x+z,x+y+z]
    for i in range(2,n):
        I = ideal(x^i, y^i, z^i, x*y*z)
        L = get_list_info(I,P)
        I1 = I.gens()
        D[I1] = L
        print((I1),L[0])
        for list in L[1]:
            print(list)
        print(L[2])
    return D

def I_ijk_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    I = (x^i,y^j,z^k) for 2 <= i<=j<=k < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,y+z,x+z,x+y+z]
    for i in range(2,n):
        for j in range(i,n):
            for k in range(j,n):
                I = ideal(x^i,y^j,z^k)
                L = get_list_info(I,P)
                I1 = I.gens()
                D[I1] = L
                print((I1),L[0])
                for list in L[1]:
                    print(list)
                print(L[2])
    return D



def J_ij_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    I = x^i,y^j-x^(j-1)*z,z^2, for 2 <= i <= j < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,y+z,x+z,x+y+z]
    for i in range(2,n):
        for j in range(i,n):
            I = ideal(x^i,y^j-x^(j-1)*z,z^2)
            L = get_list_info(I,P)
            I1 = I.gens()
            D[I1] = L
            print((I1),L[0])
            for list in L[1]:
                print(list)
            print(L[2])
    return D

def M_ij_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    I = (x^i,y^2− x*z,z^j,x^(i-1)-y,yz^(j-1)), for 2 <= i <= j < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,y+z,x+z,x+y+z]
    for i in range(2,n):
        for j in range(i,n):
            I = ideal(x^i, y^2 - x*z, z^j, x^(i-1) - y, y*z^(j-1))
            L.append(Hilbert_Function(I))
            L = get_list_info(I,P)
            I1 = I.gens()
            D[I1] = L
            print((I1),L[0])
            for list in L[1]:
                print(list)
            print(L[2])
    return D

def N_ijk_1_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    I = (x^i,y^j-x^(j-1)*z,z^2,x^(i-j+1)*y^(j-k),y^(j-k)*z), for 2 <= i  < n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,x+z,y+z,x+y+z]
    for k in range(2,n):
        for j in range(k+1,n):
            for i in range(j,n):
                if k <= i - j + 1:
                    I = ideal(x^i,y^j-x^(j-1)*z,z^2,x^(i-j+1)*y^(j-k),y^(j-k)*z)
                    L = get_list_info(I,P)
                    I1 = I.gens()
                    D[I1] = L
                    print((I1),L[0])
                    for list in L[1]:
                        print(list)
                    print(L[2])
    return D

def N_ijk_2_info(n):
    '''
    Returns a dictionary {i,j:hilbert func, jordan types, is AG} for the ideal
    I = (x^i,y^j-x^(j-1)*z,z^2,x^(i-j+1)*y^(j-k),y^(j-k)*z), for 2 < k < j <= i< n
    '''
    R.<x,y,z> = PolynomialRing(QQbar,3)
    D = {}
    P = [x,y,z,x+y,y+z,x+z,x+y+z]
    for k in range(2,n):
        for j in range(k + 1,n):
            for i in range(j,n):
                if k <= (i - j + 1):
                    I = ideal(x^i,y^j-x^(j-1)*z,z^2,x^(i-j+1)*y^(j-k),y^(j-k)*z)
                    L = get_list_info(I,P)
                    I1 = I.gens()
                    D[I1] = L
                    print((I1),L[0])
                    for list in L[1]:
                        print(list)
                    print(L[2])
    return D

def get_list_info(I,P):
    '''
    Get info for a given ideal and list of multipliers P. Returns a list with the hilbert function as the
    first entry, the jordan types as the second entry, and true or false regarding whether or not it is Gorenstein
    in the third entry.
    '''
    L = []
    J = []
    L.append(Hilbert_Function(I))
    for l in P:
        K = Jordan_Type(I,l)
        if K not in J:
            J.append(K)
    L.append(J)
    L.append(is_AG(I))
    return L


def dot_graph_for_sub_partitions(J):
    '''
    Takes in a jordan type and a hilbert function and determines if the jordan types
    could occur for that hilbert function. Returns true if so, returns false otherwise.
    ~ unfinished ~
    '''
    length = len(J)
    M = Matrix(J[0] + length)

    for j in J:
        for i in range(j):
        M[-1,i]


    #l=x case is only case
    M_1 = Matrices[0]
    for i in range(J[0]):
        M_1[-1,i] = 1
    Matrices[0] = M_1

    M_1 = Matrices[0]
    for i in range(J[1]):
        M_1[-2+i,1] = 1
    Matrices[0] = M_1

    M_2 = Matrices[1]
    for i in range(J[2]):
        M_2[-2+i,1] = 1

'''
Have theoretical checks in the beginning --> from possible_jordan_types
Find subpartitions given a jordan type
then compute smaller dot graphs for each subpartition (like the two variable case)
After that stack the graphs, shifting one over and check hilbert function by counting columns
'''

'''

def occuring_jordan_types(data,H):
    J = possible_jordan_types(H)
    for j in J:
        #check if j is in data, don't know how that is going to be formatted yet
        #if not, check using dots method if it is possible for that hilbert function using is_possible(J,H)
'''


'''
Streamline by making it functions ideals (specify by doing func for ideals defined by i, i and j, and i,j,k) but can't input an ideal
defined by variables that do not yet exist...
'''


'''
Try to write formula for N case1 and case2

H = [1,3,4^l,3,1] l > 4
Finish table on overleaf
Then... look at more families for which H = [1,3,4^l,3,1] l <= 4..?



What are the jordan types that commute with eachother?
classify those commuting pairs within individual ideals and see
all pairs that commute with eachother
How?


Write a function that takes in a hilbert and gives all jordan types that could show up
(diff from possible types, checks if the jordan types can exist for that hilbert)

Can build hilbert from jordan types, but can we build jordan types from hilbert function?
 - Maybe building possible_jordan_types then checking if they match with Hilbert func
 - Look at data (keep track of what ideal data comes from),
 and then use 8/7 notes to determine (using points method) if partition can occur
 - Use H = [1,3^k,1], don't worry about Gorenstein

With the goal of being able to tell when a possible_jordan_types cannot actually occur
'''
