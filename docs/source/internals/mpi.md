# MPI Internals

MPI support distributes selected decision diagram operations across MPI ranks.
The MPI build enables additional C++ functions and Python bindings for
distributed state creation, multiplication, measurement, and serialization.

When documenting MPI internals, keep user-facing installation steps in
{doc}`../mpi` and reserve this page for implementation constraints such as rank
ownership, communication patterns, and serialization behavior.
