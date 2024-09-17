## Simulator For Active STructures (Sim-FAST)

This is an educational code package for active structures, including tensegrity, 
origami, trusses, mechanisms, metamaterials, and many others. 
The following Intro Figure shows some working examples.

![alt text](https://github.com/zzhuyii/Sim-FAST/blob/main/Figures/Intro.png)

More specifically, the model in this package represents an active structure 
using nodal representations without rotational degrees-of-freedom. 
Such formulation are similar to the pseudo-rigid-body model coined in the 
pioneering research of Prof. Larry Howell. 
The ability to capture phenonmena like bending and twisting is enabled 
through using four-node and three-node rotational spring elements. 

While the modeling approach may not be suitable for all active structures 
(there probably exists no universal model for every problem), I believe the 
framework is particularly suitable for educational purpose. 
This formulation is good for demonstrating element formulation for
large deformation, implicit and explicit nonlinear solvers, 
and integration of simulation package. 
I am constructing this package for my own students, I hope this will
become their starting point to learn simulation for active structures. 
However, I believe such an effort is worth made open-access, which is 
why I uploaded this package here. 

## As of Sept 17th. 

I decided to ultilize the "live code" capability in MATLAB to construct this code. 
This is a good way of making educational code interactive.
It also supports showing the code and figures side by side, which is very helpful. 
The code is working for a decent number of examples, including active truss, thin origami, 
thick origami, and Bennett mechanism.  
Will continue working on the code when have time. 


## Using the code

When using the code, please add "00_SourceCode_Elements" and "00_SourceCode_Solver" 
to the path in MATLAB. For some examples, we can use "00_SourceCode_Elements_Advanced" 
for faster computation. 

Basic elements are derived using central-difference method while 
advanced elements are derived analytically. Advanced elements are all vectorized for
computation speed and they begin with "Vec_Elements". Basic elements begin with "CD_Elements"
representing that their derivation is based on central difference. 





