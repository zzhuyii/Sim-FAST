%% This function initialize the assembly system

function Initialize_Assembly(obj)

    obj.node.current_U_mat = zeros(size(obj.node.coordinates_mat));
    obj.node.current_ext_force_mat = zeros(size(obj.node.coordinates_mat));

    if isempty(obj.rot_spr_4N)
    else
        obj.rot_spr_4N.Initialize(obj.node)
    end

    obj.cst.Initialize(obj.node)

    %obj.bar.Initialize(obj.node)
    obj.rot_spr_4N_D.Initialize(obj.node)

end