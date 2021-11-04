function output=matrix3by3makesparse(atx,i,j,spinnumber)  
% originally atx contains matrices with size [3 3 Ndict], contains no spin
% info
% sparse matrix is made by S{i(k),j(k)}=value(k)
output=sparse(i',j',repmat(atx(:),[spinnumber,1]));