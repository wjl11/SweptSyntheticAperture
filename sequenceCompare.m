ref = load('P4-2_SA_4B.mat');
cmp = load('P4-2_SSA_manual');


refstart_i = 392;
cmpstart_i = cmp.SA_acq+1;

totEvent = length(ref.Event)-refstart_i;

for i = 0:totEvent
    
    ref.Event(refstart_i+i)
    disp('*******************')
    cmp.Event(cmpstart_i+i)
    keyboard

end