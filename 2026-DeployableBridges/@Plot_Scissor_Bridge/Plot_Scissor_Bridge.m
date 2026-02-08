classdef Plot_Scissor_Bridge < handle

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

        % Panel information for plotting
        panelConnection={}

    end

    methods
        % Plot the shape of the system with node number
        Plot_Shape_Node_Number(obj);

        % Plot the shape of the system with 3 node rot spr number
        Plot_Shape_RotSpr_3N_Number(obj);

        % Plot the shape of the system with 4 node rot spr number
        Plot_Shape_RotSpr_4N_Number(obj)

        % Plot the shape of the system with cst number
        Plot_Shape_CST_Number(obj)

        % Plot the number of bars
        Plot_Shape_Bar_Number(obj);

        % Plot the deformation animation
        Plot_Deformed_His(obj,Uhis)

        % Plot the deformed shape of the system
        Plot_Deformed_Shape(obj,U)

        % Plot the deformed shape of the system
        Plot_Shape_ActBar_Number(obj)

    end
end