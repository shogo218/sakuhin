function[BCCcode,trellis,puncpat]=BCCEncode(data,rate)
trellis=poly2trellis(7,[133 171]);
if strcmp(rate,'1/2') || isequal(rate,1/2)
    puncpat=ones(log2(trellis.numOutputSymbols),1);
elseif strcmp(rate,'2/3') || isequal(rate,2/3)
    puncpat=[1 1 1 0].';
elseif strcmp(rate,'3/4') || isequal(rate,3/4)
    puncpat = [1 1 1 0 0 1].';
elseif strcmp(rate,'5/6') || isequal(rate,5/6)
    puncpat = [1 1 1 0 0 1 1 0 0 1].';
end
BCCcode=convenc(data,trellis,puncpat);
return