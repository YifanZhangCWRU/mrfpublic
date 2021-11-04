function Msignal= ...
    flipmagnezationCPUParforEchoVer3rev(Atd,Btd,Ate,Bte,Atr1,Btr1,Atr2,Btr2,Atr3,Btr3, M0, dephasefactor1, dephasefactor2)
% here A is the MRF flip parameter for transversal, B is for vertical
% for example a dictionary size 200x6435 whereas 200 experiments were
% performed and 6435 T1/T2 combinations were used, you have expnumber=200
% and dictlength=6435.Totalspin number by default is 100 to start with and
% you need to verify size(A)=[dictlength,expnumber,3,3) and
% size(B)=[dictlength,expnumber,3)
% M0 refer to the first experiment m , by default 100 spins, and
% expindex=1, and for all T1 T2 combinations dictindex(:)
% Yifan Zhang 2019 October
expnumber =size(Atd,4);
spinnumber=size(M0,2);  % by default 100 spins, M0 should be 3 by 100
dictlength=size(Atd,3);
replicatefactor=spinnumber*dictlength; 
Msignal= zeros (dictlength,expnumber);
itemp1 = reshape (1:3*replicatefactor, 3,replicatefactor); 
itemp2 =repmat (itemp1,3,1);  
I= itemp2(:);
jtemp     = repmat (1:3*replicatefactor,3,1);  J=jtemp(:);
% replicate each A matrix (3x3 , from freeprecess) by this factor to make
% diagonal sparse matrices; I and J are indices to make 3x3 matrix sparse

II= (1:3*spinnumber*dictlength)';
JJtemp=repmat(1:spinnumber*dictlength, 3,1);
JJ=JJtemp(:);
% replicate each B  matrix (3x1 , from freeprecess ) by this factor to make
% diagonal sparse matrices; II and JJ are indices to make 3x1 matrix sparse

% replicate each rotation  matrix (3 x 3 x spinnumber from dephasing ) by this factor to make
% diagonal sparse matrices; III and JJJ are indices to make 3x3 matrix sparse

% example: dephase factor =2pi, then enter dephasefactor as -2*pi
phi1 = ((1:spinnumber)/spinnumber-0.5 )*(dephasefactor1);  % dephase -2pi
zrotmatrix1=double(zrotshotgun(phi1));
z1=Zrotmatrix3by3makesparse(zrotmatrix1, I,J ,spinnumber,dictlength);
phi2 = ((1:spinnumber)/spinnumber-0.5 )*(dephasefactor2/2);  % read  2pi (split 4pi into two parts)
zrotmatrix2=double(zrotshotgun(phi2));   % note : zrotmatrix is 4pi , one half of needed
z2=Zrotmatrix3by3makesparse(zrotmatrix2, I,J ,spinnumber,dictlength);

M0temp=double(repmat(M0(:), [dictlength,1]));
M0matrix=sparse(II,JJ,M0temp);
    atd=double((Atd(:,:,:,1)));      btd=double(Btd(:,:,:,1));
    atdmatrix=Amatrix3by3makesparse(atd, I, J,spinnumber);
    btdmatrix=row3vectormakesparse(btd,II,JJ,spinnumber);
    mtr3=atdmatrix*M0matrix+btdmatrix;  
tic

for expindex=1:expnumber
    disp (strcat('process MRF experiment', num2str(expindex),'/',num2str(expnumber)));
    clear mtd
    % step 1 inversion
%     atd=double((Atd(:,:,:,expindex)));      btd=double(Btd(:,:,:,expindex));
%     atdmatrix=Amatrix3by3makesparse(atd, I, J,spinnumber);
%     btdmatrix=row3vectormakesparse(btd,II,JJ,spinnumber);
%     if expindex==1
%         mtd=atdmatrix*M0matrix+btdmatrix;    
%     else 
    mtd=mtr3;
%     end
    clear atd btd atdmatrix btdmatrix mtr3
    %step 2 excitation and te waiting
    ate=double(Ate(:,:,:,expindex));        bte=double(Bte(:,:,:,expindex));
    atematrix=Amatrix3by3makesparse(ate, I, J,spinnumber);
    btematrix=row3vectormakesparse(bte,II,JJ,spinnumber);
    mte=atematrix*mtd+btematrix;
    clear ate bte atematrix btematrix
    % step3 dephase and read 
    mtr1=z1*mte;
    atr1=double(Atr1(:,:,:,expindex));        btr1=double(Btr1(:,:,:,expindex));
    atr1matrix=Amatrix3by3makesparse(atr1, I, J,spinnumber);
    btr1matrix=row3vectormakesparse(btr1,II,JJ,spinnumber);
    mtemp=atr1matrix*mtr1+btr1matrix;
    mtr2=z2*mtemp;
    atr2=double(Atr2(:,:,:,expindex));        btr2=double(Btr2(:,:,:,expindex));
    atr2matrix=Amatrix3by3makesparse(atr2, I, J,spinnumber);
    btr2matrix=row3vectormakesparse(btr2,II,JJ,spinnumber);
    mtr3=atr2matrix*mtr2+btr2matrix;
    % step4 getting magnetization signal
    mtr3x=(mtr3(1:3:end,:)); 
    tempx=reshape(diag(mtr3x),[dictlength,spinnumber]); % should be Mz of dictlength X spinnumber
    mtr3y=(mtr3(2:3:end,:));
    tempy=reshape(diag(mtr3y),[dictlength,spinnumber]); % should be Mz of dictlength X spinnumber
    Msignal(:,expindex)= ((-1).^expindex)* (mean(tempx,2)+1i.*mean(tempy,2));
    clear mtr2 mtr3x mtr3y mtemp atr1 atr1matrix atr2 atr2matrix btr1 btr1matrix btr2 btr2matrix
    %step5 continue dephasing, pass it to next experiment
    mtemp=z2*mtr3;
    atr3=double((Atr3(:,:,:,expindex)));    btr3=double(Btr3(:,:,:,expindex));
    atr3matrix=Amatrix3by3makesparse(atr3, I, J,spinnumber);
    btr3matrix=row3vectormakesparse(btr3,II,JJ,spinnumber);
    mtr3=atr3matrix*mtemp+btr3matrix;  % do NOT clear mtr3 , since it will pass to 
    % the next mrf experiment
    clear atr3 atr3matrix btr3 btr3matrix
end
    toc
    
    
   
