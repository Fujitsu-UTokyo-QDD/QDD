# Decision Diagram Internals

QDD represents quantum states and operators with decision diagrams. The central
C++ types are declared in `include/dd.h` and exposed to Python through
`qdd/bindings.cpp`.

The public documentation should explain stable concepts:

- vector and matrix edges
- node identity and sharing
- edge weights
- garbage collection thresholds
- measurement and state collapse

Implementation comments should focus on invariants and non-obvious algorithmic
choices rather than repeating each line of code.
