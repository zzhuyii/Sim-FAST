## Simulator For Active STructures (Sim-FAST)

I am hoping to create an educational code package for active structures, including tensegrity, origami, trusses, mechanisms, metamaterials, and many others. 
This package is still under construction.
Eventually, there will be a note associated with this package that will be published. 
Hoepfully I will get the opportunity of teaching it one day. 


## As of Jan 2nd. 

The code is working for a few examples, including active truss, thin origami, thick origami, mechanism, and knitting structures. 
My original thought is that using central difference methods for deriving the elements (internal 
forces and stiffness) can be very useful for an education package because it is fast. 
When teaching, we can quickly cover the derivation of elements and have time to focus on nonlinear 
solver. However, it looks like I do want to code a set of elements that use the exact
solution because tuning the step size for central difference can be tricky. Probably will start working on those as I get more time. 

Not particularly sure how much effort I can spend on this as the semester start. 
Will continue working on the code when have time.
