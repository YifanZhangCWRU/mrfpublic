function Msignal= SingleVoxelMRFSimulation(T1longtable,T2longtable,prepmode,magnetfield,parameterfile)
T1prime=T1longtable.*1000;
T2prime=T2longtable.*1000;
load (parameterfile)   % Note that all parameters are derived from parameterfile.
switch (magnetfield)
    case '9T'
        riseT=0.2;  % unit :ms
    case '7T'
        riseT=0.2181;
end
% D1=riseT; 
% D23=.151;
% TE1=1.732;
% D3=.6336;
% TE2=1.564;

switch prepmode
    case 'T1'
        M0=single(repmat([0; 0 ;-1],1,totalspins));
        [Atd,Btd]   = freeprecessCPURev(inversionTime,T1prime,T2prime,df,FA.*0);
        [Ate,Bte]   = freeprecessCPURev(te-pulseduration/2- (0.634+3*riseT+0.634),T1prime,T2prime,df,FA);
        [Atr1,Btr1] = freeprecessCPURev((0.634+riseT)*ones(Nex,1)+.1,T1prime,T2prime,df,FA.*0);
        [Atr2,Btr2] = freeprecessCPURev((0.634+0.1)*ones(Nex,1)+.1,T1prime,T2prime,df,FA.*0);
        [Atr3,Btr3] = freeprecessCPURev(tr./2,T1prime,T2prime,df,FA.*0);
        Msignal=flipmagnezationCPUParforEchoT1prep(Atd,Btd,Ate,Bte,Atr1,Btr1,Atr2,Btr2,Atr3, Btr3, M0, dephasefactor, readphasefactor);
%                 
%         M0=single(repmat([0; 0 ;-1],1,totalspins));
%         [Atd,Btd]   = freeprecessCPURev(inversionTime,T1prime,T2prime,df,FA.*0);  % inversion to excitation
%         [Ate,Bte]   = freeprecessCPURev(ones(Nex,1).*(.01+D1*3+TE1+D23+pulseduration/2),T1prime,T2prime,df,FA);   % Excitation to dephase
%         [Atr1,Btr1] = freeprecessCPURev((D3+riseT)*ones(Nex,1),T1prime,T2prime,df,FA.*0);   % dephase
%         [Atr2,Btr2] = freeprecessCPURev((D3+riseT)*ones(Nex,1),T1prime,T2prime,df,FA.*0);   % end of dephase to k space center
%         [Atr3,Btr3] = freeprecessCPURev((D3+riseT)*ones(Nex,1),T1prime,T2prime,df,FA.*0);   % second half of read
%         [Atr4,Btr4] = freeprecessCPURev((D3+riseT)*ones(Nex,1),T1prime,T2prime,df,FA.*0);   % rephase
%         [Atr5,Btr5] = freeprecessCPURev((4*riseT+D23+TE2+.014+.156+.01+pulseduration/2)*ones(Nex,1),T1prime,T2prime,df,FA.*0);   % end of rephase to next excitation
%         
%         Msignal=flipmagnezationCPUParforEchoT1prep(Atd,Btd,Ate,Bte,Atr1,Btr1,Atr2,Btr2,Atr3, Btr3, M0, dephasefactor, readphasefactor);
%         
        
    case 'T2'
        M0=single(repmat([0; 0 ;1],1,totalspins)); 
        [Atd3,Btd3]   = freeprecessCPURev(inversionTime,T1prime,T2prime,df,FA.*0);
        [Ate,Bte]   = freeprecessCPURev(te-pulseduration/2- (0.634+3*riseT+0.634),T1prime,T2prime,df,FA);
        [Atr1,Btr1] = freeprecessCPURev((0.634+riseT)*ones(Nex,1)+.1,T1prime,T2prime,df,FA.*0);
        [Atr2,Btr2] = freeprecessCPURev((0.634+0.1)*ones(Nex,1)+.1,T1prime,T2prime,df,FA.*0);
        [Atr3,Btr3] = freeprecessCPURev(tr./2,T1prime,T2prime,df,FA.*0);
        Msignal=flipmagnezationCPUParforEchoT2prep(T2prime, preptime, Atd3,Btd3,Ate,Bte...
            ,Atr1,Btr1,Atr2,Btr2,Atr3, Btr3, M0, dephasefactor, readphasefactor);
    case 'T1T2'
        M0=single(repmat([0; 0 ;-1],1,totalspins)); 
        [Atd3,Btd3]   = freeprecessCPURev(inversionTime,T1prime,T2prime,df,FA.*0);
        [Ate,Bte]   = freeprecessCPURev(te-pulseduration/2- (0.634+3*riseT+0.634),T1prime,T2prime,df,FA);
        [Atr1,Btr1] = freeprecessCPURev((0.634+riseT)*ones(Nex,1)+.1,T1prime,T2prime,df,FA.*0);
        [Atr2,Btr2] = freeprecessCPURev((0.634+0.1)*ones(Nex,1)+.1,T1prime,T2prime,df,FA.*0);
        [Atr3,Btr3] = freeprecessCPURev(tr./2,T1prime,T2prime,df,FA.*0);
        Msignal=flipmagnezationCPUParforEchoT2prep(T2prime, preptime, Atd3,Btd3,Ate,Bte...
            ,Atr1,Btr1,Atr2,Btr2,Atr3, Btr3, M0, dephasefactor, readphasefactor);
        
end



        