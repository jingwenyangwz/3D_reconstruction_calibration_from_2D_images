function [vanishing_point_3D,Q_euc,P_euc] = metricReconstruction(q_data,P_Ncam_est,Q_2cam_est)
%input:
%q_data(3,npoints,ncam)
%P_Ncam_est(3,4,ncam) projection matrices
%Q_2cam_est(4,npoints) 3D homo coordinates 

[ncoord,npoints,ncam] = size(q_data);

parallel_lines = zeros(3,6,ncam);
    % we need three pairs of parallel lines
for cam = 1:ncam
    % the first pair, the three faces of the three pairs form a corner in the
    % cube
    parallel_lines(:,1,cam) =cross(q_data(:,1,cam), q_data(:,5,cam));
    parallel_lines(:,2,cam) =cross(q_data(:,2,cam), q_data(:,6,cam));

    % the second pair
    parallel_lines(:,3,cam) =cross(q_data(:,1,cam), q_data(:,2,cam));
    parallel_lines(:,4,cam) =cross(q_data(:,5,cam), q_data(:,6,cam));

    % the third pair
    parallel_lines(:,5,cam) =cross(q_data(:,1,cam), q_data(:,4,cam));
    parallel_lines(:,6,cam) =cross(q_data(:,5,cam), q_data(:,8,cam));
end
% Build a (3,3,ncam) matrix with the position (hom. coords) of the
% vanishing points in each image
vanishing_points = zeros(3,3,ncam); 
for cam =1:ncam
    % we need three vanishing points
    vanishing_points(:,1,cam) = cross(parallel_lines(:,1,cam),parallel_lines(:,2,cam));
    vanishing_points(:,2,cam) = cross(parallel_lines(:,3,cam),parallel_lines(:,4,cam));
    vanishing_points(:,3,cam) = cross(parallel_lines(:,5,cam),parallel_lines(:,6,cam));
    
    %normalize the points by the third coordinate 
    vanishing_points(:,1,cam) =vanishing_points(:,1,cam)/vanishing_points(3,1,cam);
    vanishing_points(:,2,cam) =vanishing_points(:,2,cam)/vanishing_points(3,2,cam);
    vanishing_points(:,3,cam) =vanishing_points(:,3,cam)/vanishing_points(3,3,cam);

end
vanishing_point_3D = linear_triangulation(vanishing_points,P_Ncam_est);

Q_5points = zeros(4,5); % we need 5 points, the first three should be vanishing points, then two opposite vertices for the cube
Q_5points(:,1:3) = vanishing_point_3D; %4*3 matrix
Q_5points(:,4)= Q_2cam_est(:,1); %4*1 coords for the 1st vertice 
Q_5points(:,5)= Q_2cam_est(:,7); %4*1 coords for the 7th vertice
Heuc = calc_reference_homography(Q_5points);

% Apply that transformation to the projective points and 
% matrices. Depending on which you used to obtain the transformation
Q_euc = inv(Heuc)*Q_2cam_est;
for cam = 1:ncam
    P_euc(:,:,cam) = P_Ncam_est(:,:,cam)*Heuc;
end

[K_euc,R_euc,C_euc] = CameraMatrix2KRC(P_euc);

% Visualize results. Check both possible solutions
draw_3Dcube(Q_euc);
figure;
draw_scene(Q_euc, K_euc, R_euc, C_euc(1:3,:));
draw_3D_cube_segments(Q_euc);
figure;
S=-eye(4); S(4,4)=1;
draw_scene(S*Q_euc, K_euc, R_euc, -C_euc(1:3,:));
draw_3D_cube_segments(S*Q_euc);
