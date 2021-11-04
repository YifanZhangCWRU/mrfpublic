
abc=[6290,6295,6300,6305,6320];
toplot1=abs(Msignal(abc,:))';
plot (toplot1)
legend ('[T1=1.55sec, T2=37ms]','[T2=51ms]','[T2=70ms]','[T2=90ms]','[T2=250ms]');

