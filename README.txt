The file "model.c" writes the binary file "model.dat" which contains the accelerations and torques. This may take a while. 

"test.c" runs a simulation using "model.dat"

"rigid_bodies.c" runs a simulation from scratch, calculating accelerations and torques on the fly. It is the canonical simulation; "test.c" should produce similar results.

"common.c" contains the utility functions: linear algebra operations, file IO, and some generic rigid-body physics.

"common.h" contains all common #includes, #defines, and struct definitions

"zip.py" combines xyz files together for easy viewing of multiple trajectories. The xyz format is text-based and very easy to parse. The program VMD (Visual Molecular Dynamics) is good for viewing xyz files. 

"Makefile" makes the above via the commands make model, make test, make rigid, or just make to make all. 

