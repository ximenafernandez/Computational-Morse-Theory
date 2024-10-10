from Morse_Theory.Flow import *

# Collapsability in posets

def remove_point(X, x):
    elms=[y for y in X.list() if not y == x]
    return X.subposet(elms)

def is_beat_point(X, x):
    return len(X.upper_covers(x)) == 1  or len(X.lower_covers(x)) == 1 

def core(X): 
    for x in X:
        if is_beat_point(X,x):
            X = remove_point(X, x)
            return core(X)
    return X

def is_contractible(X):
    if X.has_top() or X.has_bottom():
        return True
    return core(X).cardinality() == 1 

def F_hat(X, x):
    elms=[y for y in X.list() if X.is_greater_than(y, x)]
    return X.subposet(elms)

def U_hat(X, x):
    elms=[y for y in X.list() if X.is_greater_than(x, y)]
    return X.subposet(elms)

def is_weak_point(X, x):
    return (is_contractible(U_hat(X, x)) or is_contractible(F_hat(X, x)))

def weak_core(X):
    for x in X:
        if is_weak_point(X, x):
            X = remove_point(X, x)
            return weak_core(X)
    return X

def is_collapsible(X): 
    '''
    Greedy algorithm, doesn't check all possible collapse routes but only chooses greedily one.
    '''
    if is_contractible(X): 
        return True
    return weak_core(X).cardinality() == 1 
    
# Vertex collapses in CW-complexes

def collapse(X, pairs = None):
    '''
    X poset (face poset of regular complex)
    Collapse pairs of cells
    '''
    if pairs is None:
        pairs = []
    elements = list(X)
    for x in X.maximal_elements():
        for y in X.lower_covers(x):
            if len(X.upper_covers(y))==1:
                pairs.append([y,x])
                elements.remove(x)
                elements.remove(y)
                return collapse(X.subposet(elements), pairs)
    return (X, pairs)

def is_collapsible_vertex(X, v):
    """
    v: vertex
    X: face poset of regular CW-complex
    """
    # Link of v
    cells = X.order_filter([(v,)])
    cells.remove((v,))
    link = X.subposet(cells)
    
    #check if link is (greedily) collapsible
    (collapsed_poset, collapse_pairs) = collapse(link)
    if len(collapsed_poset)==1: # open star collapsible, check if induced flow leads to regular complex
        pairs = collapse_pairs + [[(v,), p] for p in list(collapsed_poset)]
        
        return (True, pairs)
    
    return (False, [])
    
# Generalized Morse Core

def regular_Morse_aux(K, all_cells, critical_cells=None, sigma=None, P=None):
    """
    Recursive function to perform Morse core reduction on a simplicial complex.
    
    Parameters:
    - K: regular CW-complex to be reduced.
    - critical_cells: list of critical cells
    - sigma: acyclic matching.
    
    Returns:
    - Final reduced complex and critical cells.
    """
    X = K.face_poset()
    if critical_cells is None:
        critical_cells = []  # critical_cells
    
    if sigma is None:
        sigma = []  # matching
    
    if P is None:
        P = {}

        # Initialize flow posets P(w, z) for each pair of cells w and z as entrance path posets
        for w in X:
            for z in X:
                if X.is_less_than(w, z):
                    P[(w, z)] = entrance_path_poset(X, w, z)
                else:
                    P[(w, z)] = Poset(([], [])) 
    
    # Iterate over vertices of the complex
    for v in K.vertices():
        print('checking ', v)
        (collapsible, matching) = is_collapsible_vertex(X, v)
        if collapsible:
            print(v, 'collapsible vertex')
            print('checking regularity')
            # check regularity
            
            # Initialize flow posets P(w, z) for each pair of cells w and z as entrance path posets
            P_aux = deepcopy(P)
            
            updated_pairs = []
            regular = True

            updated_pairs = []
            regular = True  # Assume the structure is regular initially

            for (x, y) in matching:
                for w in all_cells:
                    if w != x and len(P_aux[(w, y)]) != 0:  # If P(w, y) is non-empty
                        for z in all_cells:
                            if z != y and len(P_aux[(x, z)]) != 0:  # If P(x, z) is non-empty
                                # Update P(w, z) using the 'adjoin' operation
                                P_aux[(w, z)] = adjoin(
                                    P_aux[(w, z)], (x, y), P_aux[(w, y)], P_aux[(x, z)], P_aux[(w, x)], P_aux[(y, z)]
                                )
                                updated_pairs.append((w, z))

                                # Check regularity immediately after the update
                                if len(P_aux[(w, z)]) != 0 and not is_collapsible(P_aux[(w, z)]):
                                    print(f'Non-contractible P[{w}, {z}] after update.')
                                    regular = False
                                    break 
                        if not regular:
                            break  
                if not regular:
                    break 
                    
            if regular:
                # update K and sigma
                print(v, 'regular')
                sigma += matching
                K.remove_face((v,))
                
                for p in updated_pairs:
                    P[(p[0], p[1])] = P_aux[(p[0], p[1])]

                # Recursively process the new complex
                return regular_Morse_aux(K, all_cells, critical_cells, sigma, P)
        
    # No regular collapsible vertex is found, pick a vertex and its open star as critical
    if len(K.vertices()) > 0:
        # Choose a vertex to be critical
        v = next(iter(K.vertices()))  # Pick any vertex from the complex

        print(f"No regular collapsible vertex found, removing critical vertex: {v}")

        # Save the open star of the vertex before removing it
        descending_open_star = X.order_filter([(v,)])
        critical_cells += descending_open_star

        # Remove the vertex from the complex
        K.remove_face((v,))

        # Continue iterating after removing the vertex
        return regular_Morse_aux(K, all_cells, critical_cells, sigma, P)
    
    # Return the final reduced complex, saved cells, and poset modifications
    return (K, all_cells, critical_cells, sigma, P)
    
def regular_Morse_reduction(K):
    K_red, all_cells, critical_cells, sigma, P = regular_Morse_aux(deepcopy(K), K.face_poset().list())
    return K_red, critical_cells, sigma
