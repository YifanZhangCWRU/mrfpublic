
function output=row3vectormakesparse(btx,II,JJ,spinnumber)  
% originally btx contains matrices with size [3 1 Ndict], contains no spin
% info
% sparse matrix is made by S{i(k),j(k)}=value(k)
X= repmat (btx(:),[spinnumber,1]); 
output     = sparse(II,JJ,X);