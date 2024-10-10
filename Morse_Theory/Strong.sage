def is_dominated(X, v):
    """
    Check if a vertex `v` is dominated by another vertex in the face poset `X`.

    Parameters:
    - X: Face poset of a regular CW-complex K.
    - v: A tuple representing a vertex of K (equivalently, a minimal element in the poset X).

    Returns:
    - A tuple of the form (True, dominating_vertex) if `v` is dominated by another vertex in `X`.
    - (False,) if no such dominating vertex exists.
    """
    
    # Open star of `v`
    F = X.subposet(X.order_filter([v]))
    
    # Closed star of `v`
    F_open = X.subposet(X.order_ideal(F.maximal_elements()))
    
    # Candidates for a dominating vertex
    candidates = F_open.minimal_elements()
    candidates.remove(v)  # Remove `v` from the list of candidates
    
    # Check each candidate to see if it dominates `v`
    for candidate in candidates:
        dominates = True  # Assume `candidate` dominates `v`
        
        # Check if `candidate` is less than all maximal elements of `F`
        for simplex in F.maximal_elements():
            if not F_open.is_less_than(candidate, simplex):
                dominates = False
                break  # No need to check further if one condition fails
        
        # If a dominating vertex is found, return True along with the vertex
        if dominates:
            print(f"{candidate} dominates {v}")
            return (True, candidate[0])
    
    # If no candidate dominates `v`, return False
    return (False,)


def strong_Morse_reduction(X, critical_cells=None, matching=None):
    """
    Recursive auxiliary function to perform Strong Morse Theory reduction on a regular CW-complex.

    Parameters:
    - X: Face poset of a regular CW-complex K.
    - critical_cells: List to store the open stars of removed vertices.
    - matching: Directed edges in the original face poset that need to be reversed.

    Returns:
    - A tuple containing the final reduced poset `X`, the list of `critical_cells`, and the list of `matching`.
    """
    
    if critical_cells is None:
        critical_cells = []  # Initialize the list of critical cells
    
    if matching is None:
        matching = []  # Initialize the matching

    # Base case: If there are no vertices left in the complex, return the results
    if len(X.minimal_elements()) == 0:
        return X, critical_cells, matching

    # Iterate over vertices in the complex to check for dominated vertices
    for v in X.minimal_elements():
        is_dominated_result = is_dominated(X, v)
        
        if is_dominated_result[0]:  # If `v` is dominated
            dominating_vertex = is_dominated_result[1]

            # Reverse edges involving `v` and the dominating vertex
            for simplex in X:
                if v[0] in simplex and dominating_vertex in simplex:
                    # Find the face obtained after removing `v`
                    face = tuple(x for x in simplex if x != dominating_vertex)
                    
                    # Add pair to matching
                    matching.append([face, simplex])
                    
            # Remove the dominated vertex `v` and all cells containing it
            X = X.subposet([x for x in X.list() if v[0] not in x])
            
            # Recursively process the new complex
            return strong_Morse_reduction(X, critical_cells, matching)
    
    # If no dominated vertex is found, pick a vertex and mark it as critical
    if len(X.minimal_elements()) > 0:
        # Choose any vertex to be critical
        v = next(iter(X.minimal_elements()))
        print(f"No dominated vertex found, removing critical vertex: {v}")

        # Save the open star of the vertex before removing it
        descending_open_star = X.order_filter([v])
        critical_cells.extend(descending_open_star)

        # Remove the vertex and all cells containing it from the complex
        X = X.subposet([x for x in X.list() if v[0] not in x])

        # Recursively process the new complex after removing the critical vertex
        return strong_Morse_reduction(X, critical_cells, matching)

    # Return the final reduced complex, critical cells, and matching
    return X, critical_cells, matchings


def Morse_core(X, critical_cells, matching):
    """
    Compute core to reduce a regular CW-complex.

    Parameters:
    - X: The face poset of a regular CW-complex.

    Returns:
    - The core reduction of the complex after applying Morse Theory.
    """
    
    # Get the cover relations from the original complex and adjust with the reversed edges
    cover_relations = X.cover_relations()
    for x in matching:
        cover_relations.remove([x[_sage_const_0 ], x[_sage_const_1 ]])  # Remove original directed edges
        cover_relations.append([x[_sage_const_1 ], x[_sage_const_0 ]])  # # Add reversed edges
    
    # Construct the core reduction using the modified cover relations
    core = Poset([X.list(), cover_relations], cover_relations=False).subposet(critical_cells)
    
    return core


