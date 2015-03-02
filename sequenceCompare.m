ref = load('P4-2_SA_4B.mat');
cmp = load('P4-2_SSA_manual');


refstart_i = 392;
cmpstart_i = cmp.SA_acq+1;

totEvent = length(ref.Event)-refstart_i;

for i = 0:totEvent
    disp('ref')
    ref.Event(refstart_i+i)
    ref.TX(ref.Event(refstart_i+i).tx)
    disp('*******************')
    disp('cmp')
    cmp.Event(cmpstart_i+i)
    cmp.TX(cmp.Event(cmpstart_i+i).tx)
    keyboard
end