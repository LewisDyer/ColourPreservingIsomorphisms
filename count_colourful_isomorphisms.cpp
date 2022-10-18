// GENERAL PIPELINE
// Read in the pattern graph H, and generate a nice tree decomposition, returning the decomposition and its associated bags
// Find the root arbitrarily, get the colours of vertices in H at the root bag, then get a list of all subsets of vertices in G of the same size with the same set of colours.
// For each such subset:
//      Check if there's a subgraph isomorphism with the bag at the current node. If there isn't, just return 0.
//      Then consider the type of the current node in the nice tree decomposition:
//          LEAF: Just return 1!
//          JOIN: Take the product of counts of the two children, with the same K!
//          FORGET: Get the colour of the new node added in H, then sum together counts of that child when adding a new vertex to K with the same colour.
//          INTRODUCE: Get the colour of the node removed in H, then remove the vertex in K with the same colour and recurse.
//      Use this recursion to tot up the total number of extensions at the root.
//      If you ever get a value of at least 1, terminate iterating through Ks and return true. Else if there's no more Ks to check, return False

// BOOST-SPECIFIC ISSUES
// How do we store graph colours? I'm guessing some sort of property map
// How do we check graph isomorphisms? Keeping in mind I need to store the map that defines such an isomorphism as well for forget and introduce nodes.