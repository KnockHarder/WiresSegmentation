classdef StdIP
    %STDIP Class of standard image processing functions
    %   Includes various general purpose IP routine
    %   Author: Suvadip Mukherjee (sm5vp@virginia.edu)
	%   Modified by: ChanningWong (channingwong@qq.com)
    
    properties
    end
    
    methods(Static)

        function [rgbImg, grayImg, fname, pname] = readImg2D(resize)
            [fname,pname] = uigetfile('*.*','Input Data');
            rgbImg = imread(strcat(pname,fname));
            dim = ndims(rgbImg);
            
            if( resize > 0 && resize < 1 )
                rgbImg = imresize(rgbImg,resize);
            end
            
            if(dim == 3)
                grayImg = rgb2gray(rgbImg);
            else
                grayImg = rgbImg;
            end
            grayImg = im2double(grayImg);
            
        end
        
        function h = imageOverlay(im1,im2,mag,col)
           %-- Overlay im2 over im1. Generally, im2 is binary 
%            axes(par);
           bg = im1.*(~im2); 
           R = bg+col(1)*im2; G = bg+col(2)*im2; B = bg+col(3)*im2;
           rgbim2 = cat(3,R,G,B);
           imshow(rgbim2,'InitialMagnification',mag); 
        end
       
    end
    
end

