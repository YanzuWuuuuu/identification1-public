[~, act_idx4] = min(vecnorm(bottomnode_coodinfo_sorted - [0.0;0.0;0.0], 2, 1));
[~, act_idx2] = min(vecnorm(bottomnode_coodinfo_sorted - [0.15*sqrt(2)/2;0.15*sqrt(2)/2;0], 2, 1));
[~, act_idx1] = min(vecnorm(bottomnode_coodinfo_sorted - [-0.15*sqrt(2)/2;0.15*sqrt(2)/2;0], 2, 1));
[~, act_idx3] = min(vecnorm(bottomnode_coodinfo_sorted - [0;-0.15;0], 2, 1));

[~, sen_idx1] = min(vecnorm(bottomnode_coodinfo_sorted - [-0.095;0.160;0.0], 2, 1));
[~, sen_idx2] = min(vecnorm(bottomnode_coodinfo_sorted - [0.095;0.160;0.0], 2, 1));
[~, sen_idx3] = min(vecnorm(bottomnode_coodinfo_sorted - [0.095;-0.160;0.0], 2, 1));
[~, sen_idx4] = min(vecnorm(bottomnode_coodinfo_sorted - [-0.095;-0.160;0.0], 2, 1));