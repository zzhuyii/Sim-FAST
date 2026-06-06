## Simulator For Active STructures (Sim-FAST)

This is an educational and research code package for active structures, including tensegrity, 
origami, trusses, mechanisms, metamaterials, and many others. **Please ensure running in MATLAB 2024a or newer version!**

This repo stores the MATLAB version code. **Please find the Python version using the following link:**

[https://github.com/zzhuyii/Sim-FAST-PY](https://github.com/zzhuyii/Sim-FAST-PY)

The following Intro Figure shows some working examples.

![alt text](https://github.com/zzhuyii/Sim-FAST/blob/main/03_Figures_On_Git_Webpage/Introduction.png)


More specifically, the model in this package represents an active structure 
using nodal representations without rotational degrees-of-freedom. 
Such formulation are similar to the pseudo-rigid-body model coined in the 
pioneering research by Prof. Larry Howell, the lumped mass model widely used in structural dynamics,
the early origami structure simulation developed by Prof. Resch and Prof. Cristiansen, 
and the more recent bar and hinge simulation for origami by Prof. Liu and Prof. Paulino. 
The ability to capture phenonmena like bending and twisting is enabled 
through using four-node and three-node rotational spring elements. 

While the modeling approach may not be suitable for all active structures 
(there probably exists no universal model for every problem), I believe the 
framework is particularly suitable for educational purpose and many 
interesting research problems I myself have worked on. 
This formulation is good for demonstrating element formulation for
large deformation, implicit and explicit nonlinear solvers, 
and integration of simulation package. 
The formulation can be easily adopted for many active structures rapidly
for a research project, giving people great flexibility. 
Being non-black-box also make it easy to make changes to the code for 
specific problems or connect the code with optimization packages for design. 
I will not claim that the code package is a replacement to commercial 
softwares (which is not the case), but it offers very different capabilities 
that I myself does not find in commercial softwares. 


## Youtube Tutorial Videos

Please find the list of tutorial videos here:

https://www.youtube.com/playlist?list=PLLR0SJ1WSb92-_e7nog4T6sGMNHzgWRNx

The first vdieo from this series can be found here:

[![](https://img.youtube.com/vi/CvfD8SYBZcY/0.jpg)](https://www.youtube.com/watch?v=CvfD8SYBZcY&list=PLLR0SJ1WSb92-_e7nog4T6sGMNHzgWRNx)

I tried to make my group's internal training more systematic by creating these videos. I believe that this will be a great first step for my own students to get started with coding and simulation of active and adaptive structures. Hopefully, these tutorial videos will eventually evolve into a grad level class on the topics of active structures. Since I have already created these, why not make them open-access right ;)

 

## A Living Textbook 

This package is associated with  a living textbook component (living means it is growing ;) ). 
This living textbook is based on the live script function in MATLAB, where we can combine text, equations, figures, codes, functions, and results in one place. 
This is really ideal for introducing simulation of active structures because we really need all these components side by side. 
I am constructing this package for my own students, I hope this will
become their starting point to learn simulation for active structures. 
However, I believe such an effort is worth made open-access, which is 
why I uploaded this package here.

![alt text](https://github.com/zzhuyii/Sim-FAST/blob/main/03_Figures_On_Git_Webpage/Feature.png)

**Please ensure running in MATLAB 2024a or newer version!**

In the live script ("living textbook"), functions are created as we are working through examles. Thus, deifinition of functions are merged with global codes. Running such mixture of live script is only supported in 2024a or newer. In older MATLAB versions, functions must be in the back (which is not ideal when the code is integrated into textbooks). This "living-textbook" idea is only possible with the 2024a or newer MATLAB so Hey! We are definitely the state-of-the-art!

When using the code, please add "00_SourceCode_Elements" and "00_SourceCode_Solver" 
to the path in MATLAB. For some examples, we can use "00_SourceCode_Elements_Advanced" 
for faster computation. 

Basic elements are derived using central-difference method while 
advanced elements are derived analytically. Advanced elements are all vectorized for
computation speed and they begin with "Vec_Elements". Basic elements begin with "CD_Elements"
representing that their derivation is based on central difference. 


## Reference and Acknowledgement:

I would like to acknowledge discussion and support from multiple agencies and individuals who have made this code package possible. 
I want to acknowledge helpful discussion with Dr. M. Schenk on the overall code structure and textbook idea. 
I want to acknowledge the inspiration from Dr. R. Resch and Dr. H. N. Cristiansen for their pioneering work on the simulation of origami structures (Resch and Cristiansen (1970)).
I want to acknowledge the inspiration from Dr. K. Liu and Dr. G. H. Paulino for their work on using bar-and-hinge model to simulate origami structures (Liu and Paulino (2017)).
I want to acknowledge the support of my Phd advisor Dr. E. T. Filipov for supporting the development of my previous origami simulation package. 


The triangle-to-triangle penetration prevention element uses the distance calculation method and code developed by Shellshear and Ytterlid (2014). 
The code and method from Shellshear and Ytterlid (2014) is published under the MIT liscence. 
Please ensure their work is properly cited if you are using the triangle-to-triangle penetration prevention element. 

* Resch, R., Cristiansen, H.N., 1970. Kinematic folded plate system, in: Proceedings of IASS Symposium on Folded
Plates and Prismatic Structures, IASS, Vienna, Austria.
* Liu, K., Paulino, G., 2017. Nonlinear mechanics of non-rigid origami: an efficient computational approach. Proceed-
ings of Royal Society A 473. 
* Shellshear, E., & Ytterlid, R. 2014. Fast Distance Queries for Triangles, Lines, and Points using SSE Instructions. Journal of Computer Graphics Techniques (JCGT), 3(4), 86–110. 


## Release Note:

This project is released using the CC-BY 4.0 liscence. 
You are free to share, adapt, or even commercially use the contant as long as appropriate credit is given to the creator. 
Guidline about CC-BY 4.0 liscence is found here:
https://creativecommons.org/licenses/by/4.0/deed.en

## Please cite this work as: 
Yi Zhu, Simulator for Active STructures: Sim-FAST - A Living Textbook, https://github.com/zzhuyii/Sim-FAST, Accessed at XX day XX month XX year. 


