function [E,P_euc,Q_euc] = EuclideanReconstructionUsingEssentialMatrix(q_2cam,q_est,F,K)

% q_2cam are the homogeneous coordinates of  originally projected points
% q_est are the homogenous coordinates of reprojected points 

npoints=size(q_2cam,2);
% ------------------------------------------------------------------------
% 2. Visualization of ideal projectios, noisy projections and reprojected points
% ------------------------------------------------------------------------
% for k=1:2
%      figure()  
%      hold on
%      scatter(q_2cam(1,:,k),q_2cam(2,:,k),30,[0,1,0]);        % Noisy projections. Green
%      q_est = un_homogenize_coords(q_est);
%      scatter(q_est(1,:,k),q_est(2,:,k),30,[0,0,1]);    % Reprojected points after projective calibration. Blue
%      hold off
%      title(sprintf('Image %d', k));
%      % axis([-1000, 1000, -1000, 1000]);
%      daspect([1, 1, 1]);
%      pbaspect([1, 1, 1]);
% end
% hold off

% ------------------------------------------------------------------------
% 3. Obtain the essential matrix (E) from the fundamental matrix (F) and the
% intrinsic parameter matrices (K).
% ------------------------------------------------------------------------
E = K(:,:,2).'*F*K(:,:,1);

% ------------------------------------------------------------------------
% 4. Factorize the essential matrix with the 2 possible solutions for
% R. Use the function factorize_E to obtain R_est(:,:,1) and R_est(:,:,2) and T_est.
% ------------------------------------------------------------------------
[R_est,T_est] = factorize_E(E);

% ------------------------------------------------------------------------
% Save the 4 solutions (R,t) in the structures Rcam(3,3,cam,sol), T(3,cam,sol),
% where cam indicates the camera number and sol indicates the solution number (1, 2, 3 or 4).
% ------------------------------------------------------------------------
Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);
%assign value to the first camera, identical matrix for rotation, null
%vector for translation
for i=1:4
    Rcam(:,:,1,i) = eye(3,3);
    Tcam(:,1,i) = zeros(3,1);
end

%assign to the second camera
Rcam(:,:,2,1) = R_est(:,:,1);
Tcam(:,2,1) = T_est;

Rcam(:,:,2,2) = R_est(:,:,1);
Tcam(:,2,2) = -T_est;

Rcam(:,:,2,3) = R_est(:,:,2);
Tcam(:,2,3) = T_est;

Rcam(:,:,2,4) = R_est(:,:,2);
Tcam(:,2,4) = -T_est;

% ------------------------------------------------------------------------
% 5. For each solution we obtain an Euclidean solution and we visualize it.
% ------------------------------------------------------------------------
npoints = size(q_2cam,2);
Q_euc = zeros(4,npoints,4); % Variable for recontructed points
P_euc = zeros(3,4,2,4);     % Variable for projection matrices
figNo=figure;

for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc(:,:,sol) = TriangEuc(Rcam(:,:,2,sol),Tcam(:,2,sol),K,q_est);
    
    % visualize 3D reconstruction
    figure();
    draw_scene(Q_euc(:,:,sol), K, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
    
    % Compute the projection matrices from K, Rcam, Tcam
    for i = 1:2
        R_T = [Rcam(:,:,i,sol) -Rcam(:,:,i,sol)*Tcam(:,i,sol)];
        P_euc(:,:,i,sol) = K(:,:,i)*R_T;    
    end    
    
    % Obtain the re-projected points q_rep
    for k = 1:2
        for j = 1:npoints
            q_rep(:,j,k) = P_euc(:,:,k,sol)*Q_euc(:,j,sol);
        end
    end
    
    % Visualize reprojectd points to check that all solutions correspond to
    % the projected images
    q_rep = un_homogenize_coords(q_rep);
    for k=1:2
      figure(figNo); subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
      title(sprintf('Reprojection %d, image %d', sol, k));
    end
    
 

%disp(['Resudual reprojection err Euclidean reconstruction for each solution   = ' num2str( ErrorRetroproy(q_2cam,P_euc,Q_euc(:,:,sol))/2)]);
%draw_reproj_error(q_2cam,P_euc,Q_euc(:,:,sol));





end

disp('************************************* END')

