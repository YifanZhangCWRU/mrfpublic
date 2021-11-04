function Rz=zrotshotgun(phi)
% here zrotshotgun can calculate a larger z rotation matrix 
% replicat 3x3 matrix by a factor of n, whereas n = length(phi)
% Yifan Zhang modified 2019
n=length(phi);
Rz=single(zeros(3,3,n));
for i=1:n
    Rz(:,:,i)=[cos(phi(i)) -sin(phi(i)) 0;sin(phi(i)) cos(phi(i)) 0; 0 0 1];
end