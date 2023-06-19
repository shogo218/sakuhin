function [DecodeData]=BCCDecode(code,trellis,rate,puncpat)
metric=real(code+1)/2;
qmetric=quantiz(metric,[0.001,.1,.3,.5,.7,.9,.999]);
if strcmp(rate,'1/2')||isequal(rate,1/2)
    traceback=30;
elseif strcmp(rate,'2/3')||isequal(rate,2/3)
    traceback=45;
elseif strcmp(rate,'3/4')||isequal(rate,3/4)
    traceback=60;
elseif strcmp(rate,'5/6')||isequal(rate,5/6)
    traceback=90;
end
DecodeData=vitdec(qmetric,trellis,traceback,'trunc','soft',3,puncpat);
return