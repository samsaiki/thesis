import pickle

import itertools

from sage.all import *

from sage.groups.perm_gps.permgroup_element import make_permgroup_element

def small_imaginary(complex):
    if imag(complex)<.0000001 and imag(complex) > -.00000001:
        return True
    else:
        return False

def ConvertToLetters(string):
    var('a,b,c,d,e,f,g')
    if string == 0:
        return a
    elif string == 1:
        return b
    elif string == 2:
        return c
    elif string == 3:
        return d
    elif string == 4:
        return e
    elif string == 5:
        return f
    elif string == 6:
        return g

class Tri(object):
    """Creates an abstract triangulation class from snap format triangulation gluing data"""      
    def __init__(self, TetGluings):
        self.NumberOfTet = TetGluings[0] 
        self.TetGluings = {n : TetGluings[1][n] for n in range(self.NumberOfTet)}
        self.NumberOfEdges = 0
        Edges = {}
        for n in range(self.NumberOfTet):
            """Creates a dictionary for the edges of the triangulation. """
            Edges[(n,(0,1))] = [-1,-1]
            Edges[(n,(0,2))] = [-1,-1]
            Edges[(n,(0,3))] = [-1,-1]
            Edges[(n,(1,2))] = [-1,-1]
            Edges[(n,(1,3))] = [-1,-1]
            Edges[(n,(2,3))] = [-1,-1]
        self.Edges = Edges
            


#Methods to access different levels of triangulation data
        
    def getNumberOfTet(self):
        return self.NumberOfTet
    
    def getTetGluings(self):
        return self.TetGluings

    def getEdges(self):
        return self.Edges
    
    def TetGlued(self, Tet, TetFace):
        return self.TetGluings[Tet][TetFace][0]

    def Face(self,Tet,TetFace):
        return self.TetGluings[Tet][TetFace][1]
 
    def next_face(self,tet,tet_face):
        return self.TetGluings[tet][tet_face][1][tet_face]
    
##Methods to identify edges

    
    def Sign(self,Pair): #Returns 1 if a 2-tuple is ordered, -1 otherwise
        if Pair[0] < Pair[1]:
            return 1
        else:
            return -1

    def ReverseEdge(self,Pair): #Reverses the labeling of a 2-tuple
        return (Pair[1],Pair[0])
    
    def EdgeToFaces(self,CurrentEdge): #Finds what faces to look at for the gluing data for the current edge
        Faces = [0,1,2,3]
        Faces.remove(CurrentEdge[1][0])
        Faces.remove(CurrentEdge[1][1])
        self.NextFaces = Faces

    def NextEdge(self,CurrentEdge): #Checks where the current edge glues. Orders edges, but includes a sign for orientation
        """Sets self.NextEdges to the two next possible edges in the format [Edge1,Edge2]. Sets self.NextEdgesSigns to include signs""" 
        CurrentTet = CurrentEdge[0]
        self.EdgeToFaces(CurrentEdge)
        for index in range(2):
            Edge = (self.getTetGluings()[CurrentTet][self.NextFaces[index]][1][CurrentEdge[1][0]],
             self.getTetGluings()[CurrentTet][self.NextFaces[index]][1][CurrentEdge[1][1]])
            if self.Sign(Edge) == 1:
                self.NextEdges[index] = (self.getTetGluings()[CurrentTet][self.NextFaces[index]][0],Edge)
                self.NextEdgesSigns[index] = 1
            else:
                self.NextEdges[index] = (self.getTetGluings()[CurrentTet][self.NextFaces[index]][0],self.ReverseEdge(Edge))
                self.NextEdgesSigns[index] = -1

    def LabelEdge(self,CurrentEdge,index): #Label the next edge with the current edge's variable, acounting for the orientation sign and make the next edge the current edge
        self.Edges[self.NextEdges[index]] = [self.CurrentVariable[0]*self.NextEdgesSigns[index],self.CurrentVariable[1]] #label next edge
        self.CurrentVariable = [self.CurrentVariable[0]*self.NextEdgesSigns[index],self.CurrentVariable[1]]     #
        self.CurrentEdge = self.NextEdges[index]

    def EdgeVariable(self,Edge):    #Returns the variable attached to an edge
        return self.Edges[Edge]
        
    def LabelEdges(self,CurrentEdge):
        self.NextEdge(CurrentEdge)  #Find where current edge glues
        for index in range(2):
             if self.EdgeVariable(self.NextEdges[index])[-1] != self.CurrentVariable[1]: #See if we have already glued to the next edge
                self.LabelEdge(CurrentEdge,index)      #Label Next Edge and make it the current edge    
                self.EdgesLabeled = self.EdgesLabeled + 1
                self.LabelEdges(self.CurrentEdge)      #Recurse
             else:
                pass
    
    def LabelAllEdges(self):
        self.NextEdges = [0,0]
        self.NextEdgesSigns = [1,1]
        self.NextFaces = [0,0]
        self.Variables = {}      #Initialize variables to label edges
        self.EdgesLabeled = 0
        
        index = 0   #Index for number of variables used
        
        while self.EdgesLabeled < self.getNumberOfTet():
            for EdgeLabel in self.Edges.keys():
                if self.Edges[EdgeLabel][1] not in self.Variables:
                    self.Variables[index] = [index]
                    self.CurrentVariable = [1,index]
                    self.CurrentEdge = EdgeLabel
                    self.CurrentTet = self.CurrentEdge[0]           #Change current tet to next one
                    index = index + 1
                    self.LabelEdges(self.CurrentEdge)
                else:
                    pass
                
        self.NumberOfEdges = len(self.Variables)

    def getNumberOfEdges(self):
        return self.NumberOfEdges

    def number_of_vertices(self):
        self.LabelAllEdges()
        #Create a vertex list for each tet
        vertex_list = [[0,0,0,0] for i in range(self.NumberOfTet)]
        number_of_vertices
        current_tet = 0
        current_vertex = 0
        vertex_list[current_tet][current_vertex] = 1
        next_tet = 
        next_vertex = 
        while vertex_list
    
    """Orientation methods"""

    def index_from_one(self,list):
        intermediate_list = [0]*4
        for i in range(4):
            intermediate_list[i] = list[i] + 1
        return intermediate_list

    def face_pairing_orientation(self,tet,face):
        S = SymmetricGroup(4)
        return sign(make_permgroup_element(S,self.index_from_one(self.Face(tet,face))))
    
    def orientations(self):
        orientations = [0]*self.NumberOfTet
        orientations[0] = 1
        current_tet = 0
        while orientations.count(0) != 0:
            for i in range(4):
                next_tet = self.TetGlued(current_tet,i)
                orientation = self.face_pairing_orientation(current_tet,i)            
                if orientations[next_tet] == 0:
                    orientations[next_tet] = -1*orientations[current_tet]
                    current_tet = next_tet
                    break
        return orientations

    
    """Ptolemy Variety methods"""

    def times(self,a,b):
        if a == b:
            return ConvertToLetters(str(a))+'^2'
        elif a < b:
            return ConvertToLetters(str(a))+'*'+ ConvertToLetters(str(b))
        else:
            return ConvertToLetters(str(b))+'*'+ConvertToLetters(str(b))
        
    def format_sign(self,a):
        if a == 1:
            return 1
        else:
            return -1

    def format_leading_sign(self,a):
        if a == 1:
            return 1
        else:
            return -1

    def edge_letter(self,tet,i,j):
        return self.format_leading_sign(self.Edges[(tet,(0,3))][0])*ConvertToLetters(self.Edges[(tet,(i,j))][1])
        
    def edge_sign(self,tet,i,j):
        return self.format_leading_sign(self.Edges[(tet,(0,3))][0])
        
    def PtolemyEquations(self):
        equation = []
        for tet in range(self.NumberOfTet):
            equation.append(self.format_leading_sign(self.Edges[(tet,(0,3))][0]*self.Edges[(tet,(1,2))][0]) *
            ConvertToLetters(self.Edges[(tet,(0,3))][1])*ConvertToLetters(self.Edges[(tet,(1,2))][1]) +
                            self.format_leading_sign(self.Edges[(tet,(0,1))][0]*self.Edges[(tet,(2,3))][0])
                            *ConvertToLetters(self.Edges[(tet,(0,1))][1])*ConvertToLetters(self.Edges[(tet,(2,3))][1])  ==  
            self.format_leading_sign(self.Edges[(tet,(0,2))][0]*self.Edges[(tet,(1,3))][0]) *
            ConvertToLetters(self.Edges[(tet,(0,2))][1])*ConvertToLetters(self.Edges[(tet,(1,3))][1]))
        for i in range(len(equation)):
            equation[i] = equation[i].subs(b == 1)    
        return equation
    
    """Complex volume Methods"""

    def lambda_of_triangulation(self,solution):
        block_elements=[]
        sum = [0,0]
        for tet in range(self.NumberOfTet):
            block = [0,0]
            block[0] = (self.edge_letter(tet,0,3) + self.edge_letter(tet,1,2) - self.edge_letter(tet,0,2) - self.edge_letter(tet,1,3)).subs(b == 0)
            block[0] = block[0].subs(a == log(a))
            block[0] = block[0].subs(solution[0])
            block[1] = (self.edge_letter(tet,0,1) + self.edge_letter(tet,2,3) - self.edge_letter(tet,0,2) - self.edge_letter(tet,1,3)).subs(b == 0)
            block[1] = block[1].subs(a == log(a))
            block[1] = block[1].subs(solution[0])
            block_elements.append(block)
            
        print('Block elements are')
        print(block_elements)
        return block_elements

    def numerical_flattening(self,block_element):
        """test"""
        e = block_element[0]
        f = block_element[1]
        z = exp(e)
        p = n((e-log(z))/(pi*I))
        def q(sign,z_sign):
            return (f-log(sign+z_sign*z))/(pi*I)
        if e == 0 and f == 0:
            return [0,0,0]
        if small_imaginary(n(q(1,1))) == True:
            return [z,p,real(n(q(1,1)))]
        elif small_imaginary(n(q(-1,-1))) == True:
            return [z,p,real(n(q(-1,-1)))]
        elif e == 0:
            return False
        elif small_imaginary(n(q(-1,1))) == True:
            return [z,p,real(n(q(-1,1)))]
        elif small_imaginary(n(q(1,-1))) == True:
            return [z,p,real(n(q(1,-1)))]
        else:
            print('dang')
            return False    
    
    def regulator(self,flattened_block_element): #note *I so this is complex volume
        print(flattened_block_element)
        if flattened_block_element == False:
            return False
        z=flattened_block_element[0]
        p=flattened_block_element[1]
        q=flattened_block_element[2]
        if z == 0 and p == 0 and q == 0:
            return 0
        return n((dilog(z) + 1/2*(log(z)+p*pi*I)*(log(1-z)-q*pi*I)-pi**2/6)*I,digits=10) 

    def complex_volume(self,block_elements):
        volume = 0
        orientations = self.orientations()
        for tet in range(self.NumberOfTet):
            if block_elements[tet] == False:
                return False
            block_element = block_elements[tet]
            volume = (volume + orientations[tet] *
                      self.regulator(self.numerical_flattening(block_element)))  
        return volume
       
    
    """Ptolemy co-cycle methods"""
    
    def is_simplex_cocycle(self,potential_cocycle):
        allowed_cocycles = self.allowed_cocycles()
        if allowed_cocycles.count(potential_cocycle):
            return True
        else:
            return False

    def allowed_cocycles(self):
        return [[1,1,1,1],[1,1,-1,-1],[-1,1,1,-1],[1,-1,-1,1],[-1,-1,-1,-1],[-1,-1,1,1]]

    def standard_cocycle(self):
        return [1,-1,1,-1]
    
    def coherency_check(self,cocycle):
        for tet in range(self.NumberOfTet):
            for face in range(4):
                next_tet = self.TetGlued(tet,face)
                next_face = self.next_face(tet,face)
                if cocycle[next_tet][next_face] != cocycle[tet][face]:
                    return False
        return True

    def multiply_by_neg1(self,cocycle):
        return [ x * -1 for x in cocycle]
        
    def equivalent_cocycle(self,cocycle1,cocycle2):
        if (cocycle1[0]*cocycle1[3] == cocycle2[0]*cocycle2[3] and
                cocycle1[2]*cocycle1[3] == cocycle2[2]*cocycle2[3]):
            return True

        neg_cocycle2 = self.multiply_by_neg1(cocycle2)

        product = [cocycle1[i]*cocycle2[i] for i in range(4)]
        neg_product = self.multiply_by_neg1(product)

        if product == self.standard_cocycle() or neg_product == self.standard_cocycle():
            return True
        return False
    
    def equivalent_cocycle_list(self,cocycle_list_1,cocycle_list_2):
        length = len(cocycle_list_1)
        for tet in range(length):
            if self.equivalent_cocycle(cocycle_list_1[tet],cocycle_list_2[tet]) == False:
                return False
        return True

    def is_equivalent_not_in_list(self, list_of_cocycle_lists, cocycle_list):
        for list in range(len(list_of_cocycle_lists)):
            if self.equivalent_cocycle_list(list_of_cocycle_lists[list],cocycle_list) == True:
                return False
        return True
                          
    def non_equivalent_cocycles(self,list_of_cocycle_lists):
        non_equivalent_cocycles = [list_of_cocycle_lists[0]]
        for list in range(1,len(list_of_cocycle_lists)):
            if self.is_equivalent_not_in_list(non_equivalent_cocycles,list_of_cocycle_lists[list]) == True:
                non_equivalent_cocycles.append(list_of_cocycle_lists[list])
        return non_equivalent_cocycles
     
    def cocycles(self):
        choices_of_cocycle_for_each_tet = []
        for tet in range(self.NumberOfTet):
            choices_of_cocycle_for_each_tet.append(self.allowed_cocycles())
        possible_cocycles = [p for p in itertools.product(*choices_of_cocycle_for_each_tet)]
        #possible_cocycles has format [possible_i][cocycle on tet_j][+- on face_k]
        length = len(possible_cocycles)
        allowed_cocycles = []
        num_of_tet = self.NumberOfTet
        for cocycle in range(length):
            #check if cocycles are coherent
            if self.coherency_check(possible_cocycles[cocycle]) == True:
                allowed_cocycles.append(possible_cocycles[cocycle])
        non_equivalent_cocycles = self.non_equivalent_cocycles(allowed_cocycles)
        return non_equivalent_cocycles

    def psl_ptolemy_equations(self):
        equation_list = []
        cocycles = self.cocycles()
        for list in range(len(cocycles)):
            equations = []
            sigma_list = cocycles[list]
            for tet in range(self.NumberOfTet):
                sigma = sigma_list[tet]
                (equations.append(sigma[2] * sigma[3] *
                    self.format_leading_sign(self.Edges[(tet,(0,3))][0] *
                    self.Edges[(tet,(1,2))][0]) *
                    ConvertToLetters(self.Edges[(tet,(0,3))][1]) *
                    ConvertToLetters(self.Edges[(tet,(1,2))][1]) + 
                    sigma[0] * sigma[3] *
                    self.format_leading_sign(self.Edges[(tet,(0,1))][0] *
                    self.Edges[(tet,(2,3))][0]) * ConvertToLetters(self.Edges[(tet,(0,1))][1]) *
                    ConvertToLetters(self.Edges[(tet,(2,3))][1])  ==
                    self.format_leading_sign(self.Edges[(tet,(0,2))][0] *
                    self.Edges[(tet,(1,3))][0]) *
                    ConvertToLetters(self.Edges[(tet,(0,2))][1]) *
                    ConvertToLetters(self.Edges[(tet,(1,3))][1])))
            for i in range(len(equations)):
                equations[i] = equations[i].subs(b == 1)
            equation_list.append(equations)
        return equation_list

    """Volume Methods"""

    def complex_volumes(self):
        self.LabelAllEdges()
        if self.getNumberOfEdges() == 1:
            return False
        equations = self.PtolemyEquations()
        volumes = []
        var('a,b,c')
        solutions = solve(equations, a,c,d)
        for j in range(len(solutions)):
            if (len(solutions[j]) != 0 and (str(solutions[j][0]) != 'a == 0' and
                  str(solutions[j][1]) != 'c == 0' and
                     str(solutions[j][2]) != 'd == 0')):
                sol = sols[j]
                print(equations)
                print(solutions)
                block_elements = self.lambda_of_triangulation(solutions)
                volume = self.complex_volume(block_elements)
                volumes.append(volume)
        return volumes

    def psl_volumes(self):
        volumes = []
        self.LabelAllEdges()
        if self.getNumberOfEdges() == 1:
            return False
        psl_ptolemy_equations = self.psl_ptolemy_equations()
        print(psl_ptolemy_equations)
        var('a,b,c')
        for i in range(len(psl_ptolemy_equations)):
            solutions = solve(psl_ptolemy_equations[i], a,c,d)
            print(solutions)
            for j in range(len(solutions)):
                if (len(solutions[j]) != 0 and
                    (str(solutions[j][0]) != 'a == 0' and
                      str(solutions[j][1]) != 'c == 0' and
                         str(solutions[j][2]) != 'd == 0')):
                    block_elements = self.lambda_of_triangulation(solutions[j])
                    volume = self.complex_volume(block_elements)
                    volumes.append(volume)
        return volumes
                        
    """Number of Tetrahedron"""
   
    def __str__(self):
        return "Number of tetrahedron is %s.\n %s" % (self.NumberOfTet,self.TetGluings)


class String(object):

    def __init__(self,String):
        self.UnorderedString = String
        self.NumberOfTet = int(String[0])

    def getString(self):
        return self.String
    
    def OrderString(self):
        String = self.UnorderedString[0]
        for Tet in range(self.NumberOfTet):
            index = 1 + Tet*20
            for Faces in range(4):
                String = String + self.UnorderedString[index + Faces]
                for vertices in range(4):
                    String = String + self.UnorderedString[index + Faces*4 + vertices + 4]
        return String
        
            
    def StringToTet(self,String):
        Tet = [0,0,0,0] 
        for i in range(4):       
            Tet[i] = String[0+i*5:5+i*5]
        return Tet

    def StringToFaces(self,String):
        Faces = [int(String[0:1]),String[1:5]]
        return Face

    def StringToFace(self,String):
        Face = [0,0,0,0]
        for i in range(4):
            Face[i] = int(String[i])
        return Face
    
    def Tri(self):
        """Takes in string and outputs list [number of tet, tet list]"""
        """Format of tetlist = [Face1, Face2, Face3, Face4]"""
        """Format of Face = []"""
        self.String = self.OrderString()
        self.TetList = []
        for i in range(self.NumberOfTet):
            self.TetList.append(self.String[1+i*20:21+i*20])
        Tri = [self.NumberOfTet,self.TetList]
        for Tet in range(self.NumberOfTet):
            Tri[1][Tet] = self.StringToTet(Tri[1][Tet])
            for Faces in range(4):
                Tri[1][Tet][Faces] = [int(Tri[1][Tet][Faces][0:1]),Tri[1][Tet][Faces][1:5]]
                Tri[1][Tet][Faces][1] = self.StringToFace(Tri[1][Tet][Faces][1])    
        return(Tri)        

"""File reading functions"""

def triangulation_list_to_volumes(triangulation_list):
    for string in list:
        equations = []
        num_of_tet = string[0]
        if len(string) == 20*int(num_of_tet) + 1:
            tri = String(string)
            triangulation = Tri(tri.Tri())
            
            
    

def triangulation_list_to_psl_volumes(triangulation_list):
    return(volumes)

nonzero_triangulations = []

def petronio_census_to_volumes(list_number):
    file = open("/Users/sam/Documents/Code/Ptolemy Variety/triangulations{}.txt".format(str(list_number)))
    output = open('/Users/sam/Documents/Code/Ptolemy Variety/volume{}.txt'.format(str(list_number)), 'w')

    wait = 'False'          #True if the current line shouldn't be inputed
    triangulations = []     #List of triangulations
    triangulation = ''      #String for current triangulation as it's read
    numbers = [str(number) for number in range(10)]
    manifolds = []
    index = 0
    for line in file:
        if line[0:6] == 'Angles':   #Stop reading, append triangulation to list
            wait = 'True'
            index = index + 1
            triangulations.append(triangulation)
            manifold = '\n'
            manifold = manifold + line 
        elif [wait,line[0:19]] == ['True', '-------------------']:  #Start reading new triangulation
                wait = 'False'
                manifolds.append(manifold)
                triangulation = ''
        elif wait == 'True':
            manifold = manifold + line
        else:
            for i in range(len(line)):
                if line[i] in numbers:
                    triangulation = triangulation + line[i]
    triangulation_list_to_volumes(triangulations)

    output.close()
    file.close()
    


"""Test data"""
sample_triangulation = Tri([3,[[[0,[1,2,3,0]],[0,[3,0,1,2]],[1,[2,1,0,3]],[1,[2,0,3,1]]],[[0,[2,1,0,3]],
            [0,[1,3,0,2]],[2,[2,1,0,3]],[2,[2,0,3,1]]],[[1,[2,1,0,3]],[1,[1,3,0,2]],
            [2,[1,2,3,0]],[2,[3,0,1,2]]]]])

B = [3,[[[0,[1,2,3,0]],[0,[3,0,1,2]],[1,[2,1,0,3]],[1,[2,0,3,1]]],[[0,[2,1,0,3]],
            [0,[1,3,0,2]],[2,[2,1,0,3]],[2,[2,0,3,1]]],[[1,[2,1,0,3]],[1,[1,3,0,2]],
            [2,[1,2,3,0]],[2,[3,0,1,2]]]]]

figure_eight = Tri([2,[[[1,[2,0,1,3]],[1,[0,3,1,2]],[1,[1,2,0,3]],[1,[0,2,3,1]]],
                       [[0,[2,0,1,3]],[0,[0,3,1,2]],[0,[1,2,0,3]],[0,[0,2,3,1]]] ] ])

