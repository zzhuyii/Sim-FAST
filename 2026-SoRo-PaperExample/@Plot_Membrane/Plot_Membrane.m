classdef Plot_Membrane < handle

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

        % color Number
        colorNum=1;

        % Decide if we plot the active bar (cable)
        showCable=1

    end

    methods
        % Plot the shape of the system with node number
        Plot_Shape_NodeNumber(obj);

        % Plot the shape of the system with rot spr number
        Plot_Shape_SprNumber(obj);

        % Plot the shape of the system with cst number
        Plot_Shape_CSTNumber(obj)

        % Plot the shape of the system with active bar number
        Plot_Shape_ActBar_Number(obj)

        % Plot the deformation animation
        Plot_DeformedHis(obj,Uhis)

        % Plot the deformed shape of the system
        Plot_DeformedShape(obj,U)

    end
end