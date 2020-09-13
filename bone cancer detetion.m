k= imread("2.JPG");
y=imresize(k,[256 256]);
figure,
imshow(k);
title("Enchondroma affected Image")
h = ones(5,5) / 25;
I2 = imfilter(k,h);
figure,
imshow(I2), title("Average Filtered Image")
I=I2;
disp("ENTER 1 FOR FUZZY C-Means")
if ndims(I)==3
I=rgb2gray(I);
else
I=I;
end
in=1;
I=imresize(I,[256 256]);
figure,imshow(I);title("Original image");
if in==1
[M N]=size(I);
data=double(reshape(I,[size(I,1)*size(I,2) 1]));
n=3;
[center, U, obj_fcn] =fcm(data, n);
maxU = max(U);
data = data';
wholeD = zeros(size(data));
for k = 1:n
indexk = (U(k,:) == maxU);
Ik = indexk.*data;
Ik = reshape(Ik,M,N);
result{k} = Ik;
wholeG(indexk) = k-1;
end
mean_cluster_val = zeros(3,1);
for k = 1:n
mean_cluster_val(k) = mean(center(k));
end
[mean_cluster_val,idx] = sort(mean_cluster_val);
img=result{idx(end)};
img=medfilt2(img,[2 2]);
end
BW = img;
figure,
imshow(BW)
title(['FCM clustered segmented image with ',num2str(n),'bins']);
se = strel('disk',11);
erodedBW = imerode(BW,se);
figure,
imshow(erodedBW)
title('Morphological Operations – Erosion')
se = strel('disk',4);
dilatedBW = imerode(erodedBW,se);
figure,
imshow(dilatedBW)
title('Morphological Operations – Dilation')
outt = blf(k,3,0.2,0.1);
figure,imshow(outt);
title('Bilateral Filtered Image');
I=outt;
disp('ENTER 1 FOR FUZZY C-Means')
if ndims(I)==3
I=rgb2gray(I);
else
I=I;
end
in=1;
I=imresize(I,[256 256]);
figure,imshow(I);title('Original image');
if in==1
[M N]=size(I);
data=double(reshape(I,[size(I,1)*size(I,2) 1]));
n=3;
[center, U, obj_fcn] =fcm(data, n);
maxU = max(U);
data = data';
wholeD = zeros(size(data));
for k = 1:n
indexk = (U(k,:) == maxU);
Ik = indexk.*data;
Ik = reshape(Ik,M,N);
result{k} = Ik;
wholeG(indexk) = k-1;
end
mean_cluster_val = zeros(3,1);
for k = 1:n
mean_cluster_val(k) = mean(center(k));
end
[mean_cluster_val,idx] = sort(mean_cluster_val);
img=result{idx(end)};
img=medfilt2(img,[2 2]);
end
BW1 = img;
figure,
imshow(BW1)
title(['FCM clustered segmented image with ',num2str(n),'bins']);
se = strel('disk',11);
erodedBW1 = imerode(BW1,se);
figure,
imshow(erodedBW1)
title('Morphological Operations – Erosion')
se = strel('disk',1);
dilatedBW1 = imerode(erodedBW1,se);
figure,
imshow(dilatedBW1)
title('Morphological Operations – Dilation')

function [out]=blf(im,sigd,sigr,n)
imwrite(im,'in.jpg')
figure,imshow(im); title('Input Image')
 w=(2*sigd)+1;
nim=imnoise(im,'Gaussian',n);
imwrite(nim,'nim.jpg')
figure,imshow(nim); title('Noisy Image')
[row clm]=size(im);
gw=zeros(w,w);
c=ceil(w/2); 
c=[c c];
for i=1:w    
    for j=1:w
        q=[i,j];
        gw(i,j)=norm(c-q);
    end
end
Gwd=(exp(-(gw.^2)/(2*(sigd^2))));
proci=padarray(nim,[sigd sigd],'replicate');
[row clm]=size(proci);
proci=im2double(proci);
K=sigd;
L=[-K:K];
c=K+1;
iter=length(L);
ind=1;
h = waitbar(0,'Wait...');
set(h,'Name','Bilateral Fiter Processing');
for r=(1+K):(row-K)     
    for s=(1+K):(clm-K)        
            for i=1:iter
                for j=1:iter                    
                    win(i,j)=proci((r+L(i)),(s+L(j)));
                end
            end
            I=win;
            win=win(c,c)-win;
            win=sqrt(win.^2);
            Gwi=exp(-(win.^2)/(2*(sigr^2)));    
            weights=(Gwi.*Gwd)/sum(sum(Gwi.*Gwd));
            Ii=sum(sum(weights.*I));            
            proci(r,s)=Ii;
            win=[];
    end
    waitbar(r/(row-K));
end
close(h)
proci=rpadd(proci,K); 
out=im2uint8(proci);
figure,imshow(out)
imwrite(out,'out.jpg'); title('Result Image') 
orgimg=im2double(imread('in.jpg'));
nimg=im2double(imread('nim.jpg'));
dimg=im2double(imread('out.jpg'));  
disp('PSNR'), psn=PSN(orgimg,dimg)
disp('IEF'), I=IEF(orgimg,nimg,dimg)
end

function [out]=IEF(orgimg,nimg,dimg)
orgimg =im2double(orgimg); 
nimg   =im2double(nimg);
dimg   =im2double(dimg);
out=sum(sum((nimg-orgimg).^2))/(sum(sum((dimg-orgimg).^2)));
end

function [out]=PSN(orgimg,mimg)
orgimg =im2double(orgimg);
mimg   =im2double(mimg);
Mse=sum(sum((orgimg-mimg).^2))/(numel(orgimg)); %Mse = Mean square Error
out=10*log10(1/Mse);
end

function x=rpadd(R,K)
for i=1:K
    R(1,:)=[];
    R(:,1)=[];
    [ro cl]= size(R);
    R(ro,:)=[];
    R(:,cl)=[];;
end
x=R;
end