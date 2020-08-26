clear all;
ima = imread('image1_small.jpeg');
ima = rgb2gray(ima);

[m, n] = size(ima);

image = zeros(m,n,6); 
image = uint8(image);

image(:,:,1) = ima;
image(:,:,2) = rgb2gray(imread('image2_small.jpeg'));
image(:,:,3) = rgb2gray(imread('image3_small.jpeg'));
image(:,:,4) = rgb2gray(imread('image4_small.jpeg'));
image(:,:,5) = rgb2gray(imread('image5_small.jpeg'));
image(:,:,6) = rgb2gray(imread('image6_small.jpeg'));

%Np is the number of points you want to obtain.
Np =  9;

%wid is the size of the checkerboard in millimeter.
widhtP = 160; %for the small checkerboard
%widhtP = 168; % for the big checkerboard

%display = 1 when you want to show the order for the returned point coordinates.
display = 1;

[coords, ima_pattern] = get_real_points_checkerboard_vmmc(Np, widhtP, display);
close all
xy = zeros(2,9,6);

%these are the points chosen for big checkboard:
% xy(:,:,1) = [343         748        1183         403         766        1156         451         784        1135;
%           897         908         920         521         518         511         212         197         181];
% 
% xy(:,:,2) = [441         756        1117         445         777       1160         450         799        1210;
%          891         948        1013         563         599         642         208         215         229];
%  
% xy(:,:,3) = [457         869        1247         390         859      1301         298         857        1375;
%         1043         995         953         678         642         608         190         163         137];
% 
% xy(:,:,4) = [378         872        1358         417         863        1306         451         860        1261;
%          917         914         912         458         460         457          82          77          73];
%     
% xy(:,:,5) = [438         869        1226         465         895        1250         496         920        1268;
%         1043        1029        1017         586         614         636         142         206         260];
% 
% xy(:,:,6) = [532         952        1304         453         926        1321         355         899        1342;
%         1127        1010         909         708         617         538         167         122          88];


%these are the points chosen for small checkboard:
xy(:,:,1) = [ 504         892        1201         478         929        1276         447         989        1376;
        1107         998         911         741         660         602         205         202         196];

xy(:,:,2) = [618         935        1136         504         863        1088         370         787        1040;
        1187         951         803         684         557         476         109         115         119];


xy(:,:,3) = [ 525        1061        1436         486         944        1286         459         862        1171;
        1022         959         914         475         505         523          92         166         221];


xy(:,:,4) = [  330         709        1211         436         784        1235         520         845        1247;
         932         990        1062         544         545         544         224         188         143];


xy(:,:,5) = [487         744        1078         495         757        1102         493         769        1129;
         885         954        1044         589         614         653         269         254         233];


xy(:,:,6) = [   565         887        1220         426         787        1159         258         660        1079;
        1086         995         900         813         704         592         473         343         208];


%tr_image = zeros(321,428,6);

for i = 1:6
    %xy(:,:,i) = get_user_points_vmmc(image(:,:,i));

    homo{i} = homography_solve_vmmc(coords',xy(:,:,i));
    t_H = maketform('projective', homo{i}');
    
    [homo{i}, rperr] = homography_refine_vmmc(coords',xy(:,:,i), homo{i});
    t_H_refined = maketform('projective', homo{i}');
    
    tr_image(:,:,i) = imtransform(ima_pattern,t_H_refined,'XData',[1 size(image(:,:,i),2)], 'YData',[1 size(image(:,:,i),1)]);
    
    subplot(121), imshow(image(:,:,i));    
    subplot(122), imshow(tr_image(:,:,i));
    pause;
    
end
%calculate the internal matrix 
A = internal_parameters_solve_vmmc(homo);

% external matrix calculation
[R,T] = external_parameters_solve_vmmc(A, homo);

%calculate the ange between the image axis
ange = radtodeg(acos (-1 * A(1,2)/A(1,1)));

