# (C) (D)ifferentiable (O)perations (L)ibrary

CDOL is a library for representing arbitrary computation as a directed
hypergraph and then executing forward and adjoint through this graph.
Conceptually, each vertex in the hypergraph corresponds to an argument to some
operation, and each hyperedge corresponds to a function taking an arbitrary
number of inputs and returning and arbitrary number of outputs.

## Features
- [x] Basic API for vertices, hyperedges, and graph creation
- [x] Source/sink detection
- [x] Forward/Adjoint execution through arbitrary hypergraphs
- [ ] Representation of graph in vertex API (i.e. recursive graph execution)
- [ ] MPI support through graph partitioning
- [ ] Python interface
