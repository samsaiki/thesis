import csv
import itertools
import collections

Linear_Dependence = collections.namedtuple('Linear_Dependence',
                                           ['triangulation',
                                            'volume','lindep'])
Manifolds = collections.namedtuple('Manifolds',['Names','Volumes'])
Manifold = collections.namedtuple('Manifold',
                                  ['triangulation','volumes','number_fields'])

def load_nonzero_volumes(list_of_nonzero_volumes):
    file = open(list_of_nonzero_volumes)
    triangulations = pickle.load(file)
    file.close()
    return(triangulations)
    

def manifolds_with_invariant_trace_field(polynomial):
    manifold_names = []
    volumes_list = []
    with open('all_volumes.csv')as csvfile:
        filereader = csv.DictReader(csvfile, delimiter=";")
        for manifold in filereader: 
            if polynomial == manifold['InvariantTraceField']:
                volume = manifold['Volume']
                if sage_eval(volume).n(digits=61) in volumes_list:
                    pass
                else:
                    volumes_list.append(sage_eval(volume).n(digits=61))
                    manifold_names.append(manifold['Name'])
    manifolds = Manifolds(manifold_names,volumes_list)
    output = open('manifolds_polynomials', 'w')
    pickle.dump(manifolds,output)
    output.close()
    return manifolds                                       

def is_linear_combination(volume_polynomial):
    volume = volume_polynomial[0]
    polynomial = volume_polynomial[1]
    manifolds_with_same_poly = manifolds_with_invariant_trace_field(polynomial)
    lindep_list = [volume] + manifolds_with_same_poly.Volumes
    pari_list = pari(lindep_list)
    return [manifolds_with_same_poly[0],pari_list.lindep()]

def manifold_to_nonzero_volumes(triangulation):
    regina_manifold = tri_list_to_regina(element)
    number_fields = NTriangulationForPtolemy(regina_manifold).ptolemy_variety(2,'all').compute_solutions().number_field()
    solutions = NTriangulationForPtolemy(man1).ptolemy_variety(2,'all').compute_solutions()
    volumes = solutions.complex_volume_numerical(drop_negative_vols=True)
    nonzero_volumes = []
    nonzero_fields = []
    for index,obstruction_class in enumerate(Vols):
        if len(obstruction_class) == 0:
            pass
        elif len(obstruction_class) == 1:
            real_volume = obstruction_class[0][0].real()
            if round(real_volume,61) == 0.0:
                pass
            else:
                nonzero_volumes.append(real_volume.n(digits=61))
                nonzero_fields.append(number_field[index][0])
        else:
            print('more than one volume')
    nonzero_volumes_fields = Manifold(triangulation,nonzero_volumes,nonzero_fields)
    return nonzero_volumes_fields
    
def manifold_list_to_lindep(triangulation_list):
    checked_volumes = []
    list_of_dependences = []
    for manifold in triangulation_list:
        nonzero_volumes = manifold_to_nonzero_volumes(manifold)
        volume = nonzero_volumes.volumes[0].n(digits=61)
        if volume in checked_volumes:
            pass
        else:
            number_field = nonzero_volumes.number_fields[0]
            linear_dependence = is_linear_combination([volume,number_field])
            manifold_linear_dependence = Linear_Dependence(manifold,
                                                           volume,
                                                           linear_dependence)
            list_of_dependences.append(manifold_linear_dependence)
            checked_volumes.append(volume)
    output = open('linear_dependences', 'w') 
    pickle.dump(list_of_dependences,output)
    output.close()
    return list_of_dependences

def load_lindep(lindep_list):
    file = open('lindep_list')
    lindep_list = pickle.load(file)
    file.close()
    return lindep_list    
    
