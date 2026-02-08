%% This function initialize the assembly system

function Initialize_Assembly(obj)

    obj.node.current_U_mat = zeros(size(obj.node.coordinates_mat));
    obj.node.current_ext_force_mat = zeros(size(obj.node.coordinates_mat));

    obj.bar.Initialize(obj.node);
    obj.actBar.Initialize(obj.node);
    obj.rot_spr_4N.Initialize(obj.node);
    obj.cst.Initialize(obj.node);

end