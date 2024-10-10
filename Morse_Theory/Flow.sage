from Morse_Theory.Flow import *

def entrance_path_poset(X, w, z):
    """
    Construct the entrance path poset `Ent[X](w, z)`.
    
    The entrance path poset is the poset of chains in `X` between `w` and `z`, ordered by inclusion.

    Parameters:
    - X: A face poset of a regular CW-complex K.
    - w: The starting element in the chains X.
    - z: The ending element in the chains X.

    Returns:
    - entrance_path_poset: The entrance path poset of chains between `w` and `z`.
    """
    
    # Extract the interval subposet between `w` and `z`
    I = X.subposet(X.interval(w, z))

    # Generate chains within the interval, filtering those that start at `w` and end at `z`
    chains = [tuple(chain) for chain in I.chains() if len(chain) >= 2 and chain[0] == w and chain[-1] == z]

    # Construct the entrance path poset using the filtered chains, ordered by inclusion
    entrance_path_poset = Poset((chains, lambda c1, c2: Set(c1).issubset(Set(c2))))
    
    return entrance_path_poset


def cartesian_path_posets(P, Q):
    """
    Construct the Cartesian product of two path posets `P` and `Q`.
    
    The elements of the product poset are formed by concatenating elements from `P` and `Q`. 
    Relations are established based on the covers of elements in `P` and `Q`.

    Parameters:
    - P,Q: path posets.

    Returns:
    - The Cartesian product poset.
    """
    # Form Cartesian product elements by concatenating elements from `P` and `Q`
    elements = [p + q for p in list(P) for q in list(Q)]

    # Establish relations based on upper covers of elements in `P` and `Q`
    relations = [[p + q, pp + qq] 
                 for p in list(P) for q in list(Q) 
                 for pp in P.upper_covers(p) + [p] 
                 for qq in Q.upper_covers(q) + [q] 
                 if p + q != pp + qq]
    
    return Poset([elements, relations])


def adjoin(P_wz, pair, P_wy, P_xz, P_wx, P_yz):
    """
    Augment the poset `P(w, z)` by including paths that arise from inverting the pair `(x > y)`, and identify elements.

    Parameters:
    - P_wz: The original poset `P(w, z)`.
    - pair: A matching pair `(x, y)` to be inverted.
    - P_wy: Poset `P(w, y)`.
    - P_xz: Poset `P(x, z)`.
    - P_wx: Poset `P(w, x)`.
    - P_yz: Poset `P(y, z)`.

    Returns:
    - updated_P_wz: The new poset `P(w, z)` with the augmented paths and identified elements.
    """
    (x, y) = pair

    # Step 1: Create the Cartesian product poset N = P(w, y) × P(x, z)
    N = cartesian_path_posets(P_wy, P_xz)
    
    # Step 2: Initialize variables for paths and relations
    N_non_identified = list(N)
    new_relations = []

    # Step 3: Identify paths in `P(w, z)` resulting from (x > y) inversion

    # Iterate through paths in `P(w, x)` and `P(x, z)` to create identified paths
    for gamma in P_wx:
        for delta in P_xz:
            # Create paths and relations involving (x > y) inversion
            path1 = gamma + delta[1:]  # Path in `P(w, z)`
            path2 = gamma + (y,) + delta  # Path in `P(w, y) × P(x, z)`

            # Establish relations using the upper and lower covers of path2
            new_relations += [[path1, path] for path in N.upper_covers(path2)] + [[path, path1] for path in N.lower_covers(path2)]

            # Remove path2 from non-identified elements if it's in N
            if path2 in N_non_identified:
                N_non_identified.remove(path2)
    
    # Repeat the process for paths in `P(w, y)` and `P(y, z)`
    for gamma in P_wy:
        for delta in P_yz:
            path1 = gamma + delta[1:]  # Path in `P(w, z)`
            path2 = gamma + (x,) + delta  # Path in `P(w, y) × P(x, z)`

            new_relations += [[path1, path] for path in N.upper_covers(path2)] + [[path, path1] for path in N.lower_covers(path2)]
            if path2 in N_non_identified:
                N_non_identified.remove(path2)

    # Step 4: Create the new poset by combining existing elements and relations with new ones
    elements = list(P_wz) + N_non_identified
    relations = P_wz.cover_relations() + new_relations + N.subposet(N_non_identified).cover_relations()

    # Construct the updated poset with the new elements and relations
    updated_P_wz = Poset((elements, relations))

    return updated_P_wz


def flow_poset(X, Sigma):
    """
    Compute the flow poset FloΣ[X] for a CW complex `X` given a Morse matching `Sigma`.

    Parameters:
    - X: Face poset of a regular CW-complex.
    - Sigma: Morse matching for the complex `X`.

    Returns:
    - A dictionary mapping pairs `(w, z)` to their corresponding flow poset `P(w, z)`.
    """
    # Dictionary to store flow posets for each pair `(w, z)`
    P = {}

    # Step 1: Initialize flow posets `P(w, z)` for each pair of cells `w` and `z` as entrance path posets
    for w in X:
        for z in X:
            if X.is_less_than(w, z):
                P[(w, z)] = entrance_path_poset(X, w, z)
            else:
                P[(w, z)] = Poset(([], []))  # Empty poset for pairs that don't satisfy `w < z`

    # Step 2: Apply the flow algorithm based on the Morse matching `Sigma`
    for (x, y) in Sigma:
        print(f"Processing matching pair: ({x}, {y})")
        for w in X.list():
            if w != x and len(P[(w, y)]) != 0:  # If `P(w, y)` is non-empty
                for z in X.list():
                    if z != y and len(P[(x, z)]) != 0:  # If `P(x, z)` is non-empty
                        # Use `adjoin` operation to update `P(w, z)`
                        P[(w, z)] = adjoin(P[(w, z)], (x, y), P[(w, y)], P[(x, z)], P[(w, x)], P[(y, z)])
    
    return P
