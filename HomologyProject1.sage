

import re
import ast
# Importing regex for cleaning file

def cousinlist():
#Gives the list of cousins in file K7family.py

  filename = "K7family.py"
  file = open(filename, "r") #opens the file

  data = file.read() #reads the file

  #print("Data:", data)

  Cousins = data.split("Cousin") #Creates a list with each element split by string "Cousin"

  return Cousins

def MMIK():
#Gives the 14MMIK. 9, 14, 16, 17, 19, and 20, are not IK. The remaining 14 are MMIK.
  list = cousinlist() #calls the cousinlist function
  notIK=[0,9-1,14-2,16-3,17-4,19-5,20-6] #list of element indices that are not MMIK
  for i in notIK:
    list.pop(i) #takes out the list element at index of list
  return list #gives a list of MMIK graphs

def converttolist(liststring):
  #replaces elements in the list with string { to [
  newlist=[]
  for l in liststring:
    string=l.replace("{","[")
    newlist.append(string)
  return newlist


def converttolist2(liststring):
  #replaces elements in the list with string } to ]
  newlist=[]
  for l in liststring:
    string=l.replace("}","]")
    newlist.append(string)
  return newlist

def convertString2List(str):
  #converts the string of cousins to the string with only the string of edges
  x = re.findall(r'\[\d+\,\d+\]', str)
  for i in range(len(x)):
    x[i] = ast.literal_eval(x[i])
  return x

def convertAll2List(lst):
  #converts the list of cousins to the list with only the list of edges
  for i in range(len(lst)):
    lst[i] = convertString2List(lst[i])
  return lst

#What to call in sage to be able to proceed with the functions below
#x=cousinlist()
#y=converttolist(x)
#z=converttolist2(y)
#l=convertAll2List(z)
#l.pop(0)

def computerhomology(lst):
#calculates the homology of the clique complex given a list for a specific graph
  dictionary={}
  for i in range(len(lst)):
    gra=Graph(lst[i])
    c=gra.clique_complex()
    homology=c.homology()
    dictionary["Cousin", i+1]=homology
  return dictionary

def simplicialhomology(lst):
#calculates the homology of the simplicial complex given a list for a specific graph
  dictionary={}
  for i in range(len(lst)):
    gra=SimplicialComplex(lst[i])
    homology=gra.homology()
    dictionary["Cousin",i+1]=homology
  return dictionary

def computefacets(lst):
#computes the facets of a graph given a list for the specific graph
  dictionary={}
  for i in range(len(lst)):
    gra=Graph(lst[i])
    c=gra.clique_complex()
    face=c.facets()
    dictionary["Cousin",i+1]=face
  return dictionary

def computerdegree(lst):
#computes the degrees of each of the vertices of the graph
  dictionary={}
  #dict2={}
  for i in range(len(lst)):
    gra=Graph(lst[i])
    foo=max(lst[i])
    print(foo)
    foo2=foo[1]
    print(foo2)
    #print(type(foo2))
    for j in range(foo2):
      print(j+1)
      degr=gra.degree(j+1)
      print("d",degr)
      dictionary[j]=degr
      #dict2[i]=dictionary
  return dictionary

def MatchingComplex(g,degree):
  graph=Graph(len(g.edges())+1)
  graph.delete_vertex(0)

def IsMatching(G,U,n):
#computes whether a subgraph U in G is a part of the Matching Complex of degree n
  if U.is_subgraph(G):
    large=max(U.edges(labels=False))[1]
    for i in range(large):
      if U.degree(i+1)<=n:
        return True
      else:
        return False
  else:
    return False

def dictdeg(list):
  dict={}
  graph=Graph(list)
  for num in graph.vertices():
    dict[num]=graph.degree(num)
  return dict
