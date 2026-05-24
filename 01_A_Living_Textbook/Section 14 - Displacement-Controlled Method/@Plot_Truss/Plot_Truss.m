
classdef Plot_Truss < handle

    properties
        % Assembly of structure
        assembly

        % Control the range for plotting
        viewAngle1=45;
        viewAngle2=45;
        displayRange=1;
        displayRangeRatio=0.2;

        % Figure size and location control
        width=800;
        height=600;
        x0=0;
        y0=0;

        % hold time for gif
        holdTime=0.01;

        % the name for animation file
        fileName='animation.gif'        

        % number of active truss (will be plotted with different color)
        activeTrussNum

    end

    methods
        % Plot the shape of the system with node number
        Plot_Shape_Node_Number(obj);

        % Plot the shape of the system with bar number
        Plot_Shape_Bar_Number(obj);

        % Plot the shape of the system with spring number
        Plot_Shape_Spr_Number(obj);

        % Plot the deformation animation
        Plot_Deformed_His(obj,Uhis)

        % Plot the deformed shape of the system
        Plot_Deformed_Shape(obj,U)

        % Plot the force of the bar
        Plot_Bar_Force(obj,F)

    end
end