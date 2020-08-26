function DoG  = extractDoGGivenSigmas(SS,sigmas)

DoG = zeros(size(SS{1},1),size(SS{1},2),numel(SS));

for j=1:(numel(SS)-1)
   
    %dggxx  = ndgauss([size(SS{j},1) size(SS{j},2)],[sigmas(j),sigmas(j)],'der',[2 0]);
    %Lxx    = imfilter(SS{j},dggxx,'symmetric','same');
    
    %dggyy  = ndgauss([size(SS{j},1) size(SS{j},2)],[sigmas(j),sigmas(j)],'der',[0 2]);
    %Lyy    = imfilter(SS{j},dggyy,'symmetric','same');

    %normf   = sigmas(j).^(2)/(sigmas(j).^(2)-sigmas(j-1).^(2));% t
    delta = sigmas(j+1)- sigmas(j);
    normf   = (sigmas(j+1)/delta).^2;% t

    DoG(:,:,j) = normf.*(SS{j+1} - SS{j});
    
end