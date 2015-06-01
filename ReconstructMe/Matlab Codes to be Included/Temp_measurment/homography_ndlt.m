function [H]=homography_ndlt(X1,X2)

% compute the centroid of X1 and X2
centrX1= [mean(X1(1:2,:),2);0];
centrX2= [mean(X2(1:2,:),2);0];

% transform the origin of points X1 and X2 to their centroid 
for i=1:length(X1);
X1_c(:,i)= X1(:,i) - centrX1;
X2_c(:,i)= X2(:,i) - centrX2;
L1(:,i)= norm(X1_c(1:2,i),2);
L2(:,i)= norm(X2_c(1:2,i),2);
end

mean_L1= mean(L1);
mean_L2= mean(L2);
s1= sqrt(2)/mean_L1;
s2= sqrt(2)/mean_L2;
%% transform the X1 and X2 to the normalized X1_n and X2_n (average length= sqrt(2)) 
for i=1:length(X1_c);
X1_n(:,i)= [s1*X1_c(1:2,i);1];
X2_n(:,i)= [s2*X2_c(1:2,i);1];
end

%% obtain the transformation T1 and T2
s1I= eye(2)*s1;
t1= [-s1*centrX1(1);-s1*centrX1(2);1];
T1 = cat(2,[s1I;0 0 ],t1);

s2I= eye(2)*s2;
t2= [-s2*centrX2(1);-s2*centrX2(2);1];
T2 = cat(2,[s2I;0 0 ],t2);

A=[];
%% Applying DLT algorithm to find H_n
for i = 1: length(X1);  
xi=X1_n(1,i);yi=X1_n(2,i);zi=X1_n(3,i);    
xi_p=X2_n(1,i);yi_p=X2_n(2,i);zi_p=X2_n(3,i); 

Ai=[0 0 0 -zi_p*xi -zi_p*yi -zi_p*zi yi_p*xi yi_p*yi yi_p*zi;
    zi_p*xi zi_p*yi zi_p*zi 0 0 0 -xi_p*xi -xi_p*yi -xi_p*zi];
A= cat(1,A,Ai);
end

%get the SVD of matrix A
[U,S,V]  =svd(A);
% h is the last row of V
h_n= V(:,9);

% Assemble the homography matrix H_n for normalized points X1_n & X2_n
H_n= [h_n(1:3,:)';h_n(4:6,:)';h_n(7:9,:)'];

% Obtains the Homography H by the formula : H= (T2)^-1 * H_n * T1
H= inv(T2)*H_n*T1;


 

