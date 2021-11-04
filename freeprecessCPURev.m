function [Afp,Bfp]=freeprecessCPURev(T,T1,T2,df,FA)
% ---------- Yifan revised 4/25/2018 CPU only
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in ms, off-resonance in Hz.
% here T ,df and FA size are dependent on MRI experiment number
% T is allowed freeprecess time, df is offresonance (Hz) and FA is flip
% angle (radian)
% while T1 and T2 are dependent on dictionary size
% Note : FA flipmatrix is for MRF specific application of precessing
% in case there's variable excitation flip angle to calculate Ate
% however if no excitation is excuted before freeprecessing, you can set 
% FA=0, that willl make rflip matrix a 3x3 identify matrix
% 
% tic
Nex=length(T);
Ndict=length(T1);
zrotmatrix=zeros(3,3,Nex);
alternate=ones(size(FA));    % Yifan added alternatve index for -1 1 -1....

% assuming size of T (related to exp number) << size (T1),aka dictionary
% size, reshape T to do matrix multiplication to speedup calculation
T0=T';
E1=exp(-(1./T1)*T0);    % equivalent version of exp(T/T1), manipulated matrix
E2=exp(-(1./T2)*T0);

% tic
for expIndex=1:Nex
%     arrange FA matrix according to indices, 
%         Even number                     odd number
%         [cosFA  0 sinFA                 [cosFA 0 -sinFA
%             0   1    0                   0     1   0   
%         -sinFA  0 cosFA]                sinFA  0   cosFA]       
    alternate(expIndex)=(-1)^expIndex; 
    zrotmatrix(:,:,expIndex)=zrot(2*pi.*df(expIndex).*T(expIndex)/1000); 
    %zrotmatrix is diagonal 3x3 matrix when off resonance df is zero
    % if df is not zero then the turn zrotmatrix(1,2) and zrotmatrix(2,1)
    % are nonzero
end
% toc

% define each zrotmatrix elements 3x3 matrix as follows
% [z1 z2 0
%  z3 z4 0
%  0   0 z5]

E2=E2';E1=E1';
E2vector=E2(:);
E1vector=E1(:);
% define cosFA 
FAflipvector1=cos(FA);
FAflipvector2=sin(FA).*alternate;
zvector1= squeeze(zrotmatrix(1,1,:));
zvector2= squeeze(zrotmatrix(1,2,:));
zvector3= squeeze(zrotmatrix(2,1,:));
zvector4= squeeze(zrotmatrix(2,2,:));
zvector5= squeeze(zrotmatrix(3,3,:));

FA1=repmat(FAflipvector1,[1 Ndict]); FA1=FA1(:);
FA2=repmat(FAflipvector2,[1 Ndict]); FA2=FA2(:);
z1=repmat(zvector1,[1 Ndict]); z1=z1(:);
z2=repmat(zvector2,[1 Ndict]); z2=z2(:);
z3=repmat(zvector3,[1 Ndict]); z3=z3(:);
z4=repmat(zvector4,[1 Ndict]); z4=z4(:);
z5=repmat(zvector5,[1 Ndict]); z5=z5(:);

% define AFP matrix elements like this
%    [ A1 A2 A3
%      A4 A5 A6
%      A7 A8 A9]    

A1=E2vector.*z1.*FA1;
A2=E2vector.*z2;
A3=E2vector.*z1.*FA2;
A4=E2vector.*z3.*FA1;
A5=E2vector.*z4;
A6=E2vector.*z3.*FA2;
A7=-E1vector.*z5.*FA2;
A8=zeros(size(E1vector));
A9=E1vector.*z5.*FA1;

% Note that FA2 is the sin with alternatinv sign ; FA1 is the cos

B1=1-E1vector;
Afp=single([A1'; A2'; A3'; A4'; A5';A6';A7';A8';A9']);
Bfp=single([B1']);

% B=[0 
%    0
%    1-E1vector]

Afp=reshape (Afp, [3,3,Nex, Ndict]); 
Afp=permute (Afp, [2 1 4 3]);  
% transpose Afp to comform matlab indexing ; also move the expeirmentnumber
% to the last index. [3 x 3 x Ndict X Nex ]
temp=reshape (Bfp, [Nex, Ndict]); 
Bfp= zeros (3,1,Ndict,Nex);
Bfp (3,1,:,:)=temp';
% toc
