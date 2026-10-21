function [num,belta,info] = SRIVC1(f,frf,Weight,belta0,noni,e,Max_iter,isubsystem)
% Joint rigid/first-flexible rational fit, with two fixed poles at zero.
assert(isubsystem==1,'This stage fits one combined rigid/flexible subsystem.');
assert(size(frf,1)==noni,'FRF channel count mismatch.');
[num,belta,info]=riv_frequency_fit(f,frf,Weight,belta0,2,2,e,Max_iter,[]);
end
