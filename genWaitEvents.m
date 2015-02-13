function [Event,n] = genWaitEvents(Event,n,scIdx,nWait)
% generate nWait number of noop events given the appropriate sequence
% control number

for i = 1:nWait
    Event(n).info = 'Wait'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;       
    Event(n).recon = 0;      
    Event(n).process = 0;
    Event(n).seqControl = scIdx; 
    n = n+1;
end