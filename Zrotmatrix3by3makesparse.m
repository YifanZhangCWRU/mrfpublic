function output=Zrotmatrix3by3makesparse(ztx,i,j,spinnumber,dictlength)  
% originally atx contains matrices with size [3 3 Ndict], contains no spin
% info
% sparse matrix is made by S{i(k),j(k)}=value(k)
% Yifan Zhang 2019 Oct
ztx=repmat(reshape(ztx,[9 spinnumber]), [dictlength,1]);              
output=sparse(i',j',ztx(:));