classdef Plot_Origami < handle

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

        % Animation file name
        fileName='Animation.gif'

        % panel information for plotting
        panelConnection={}

    end

    methods
        % Plot the shape of the system with node number
        Plot_Shape_Node_Number(obj);

        % Plot the shape of the system with bar number
        Plot_Shape_Bar_Number(obj);

        % Plot the deformation animation
        Plot_Deformed_His(obj,Uhis)

        % Plot the deformed shape of the system
        Plot_Deformed_Shape(obj,U)

    end
end