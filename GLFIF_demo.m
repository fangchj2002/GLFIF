%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Fuzzy region-based active contour driven by global and local
%  fitting energy for image segmentation" submitting to Applied Soft
%  Computing 
% Jiangxiong Fang(fangchj2002@163.com)
% East China University of Technology
% 6th, August, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
addpath 'images'
ImgID = 2;
Img = imread([num2str(ImgID),'.bmp']);

tic;

%setting the initial level set function 'u':
[M,N,L] = size(Img);
u = zeros(M,N);

u(:,:) = 0.3;
u(40:60,40:60) = 0.7;

switch ImgID
    case 1
        Img_gray = Img;
        rad = 3;
        sigma = 4;
        iterNum = 100;
        lambda1 =1.1;
        lambda2 = 1;
        alpha1 = .01;
        alpha2 = .01;
    case 2
        Img_gray = rgb2gray(Img);
        iterNum = 100;
        rad = 3;
        sigma = 3;
        lambda1 = 1;
        lambda2 = 1;
        alpha1 = 1;
        alpha2 = 1;
    case 3
        Img_gray = Img;
        iterNum = 100;
        rad = 3;
        sigma = 3;
        lambda1 = 1.8;
        lambda2 = 1;
        alpha1 = 0.3;
        alpha2 = 0.3;
    case 4
        Img_gray = Img;
        iterNum = 100;
        rad = 3;
        sigma = 3;
        lambda1 = 1.5;
        lambda2 = 1;
        alpha1 = .1;
        alpha2 = .1;
end

[Ix,Iy] = gradient(double(Img_gray));
f = abs(Ix)+abs(Iy);
g = 1./(1+f);  % edge indicator function.

figure;subplot(2,2,1);imshow(Img);hold on;%axis off,axis equal
title('Initial contour');
[c,h] = contour(u-0.5,[0 0],'r');

subplot(2,2,2);

% sigma = 5;

Ksigma = fspecial('gaussian',sigma,1.5); % Caussian kernel    

pause(0.1);
% start level set evolution
for n=1:iterNum
    [u,deltaF] = GLFIF(double(Img_gray),double(Img_gray),u,Ksigma,lambda1,lambda2,alpha1,alpha2,g);     
    if mod(n,5)==0
        pause(0.1);
        imshow(Img, []);hold on;axis off,axis equal
        [c,h] = contour(u-0.5,[0 0],'r');
        iterNum=[num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
end
subplot(2,2,3);
seg = ((u-0.5)>0);
imshow(seg);
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);




