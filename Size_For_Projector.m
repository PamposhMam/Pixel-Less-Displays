% Step 1: Read and Resize the Image
target = imread("Final_Hologram_a1_15.bmp");
%target = imresize(target, [1280, 1080]);  % Resize correctly
%target = im2bw(target);  % Convert to binary (Assuming edge image)
target=imresize(target, [1024, 1280], "bilinear");
imwrite(target, 'STarHologram_For_Projector.bmp')
imshow(target);