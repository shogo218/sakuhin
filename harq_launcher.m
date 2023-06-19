%%IR-HARQを実現するための主に符号化処理を行う関数
function [code,trellis,puncpat] = harq_launcher(data,cfgHARQ,mode,...
    psdu)

rate = cfgHARQ.rate;

switch mode
    %初送時の処理、符号化
    case 1
        %特定の符号化率の時にバッファを確保
        if strcmp(rate,'5/6')||isequal(rate,5/6)
            blank = 5;
            blankSize = mod(length(data),5);
        elseif strcmp(rate,'3/4')||isequal(rate,3/4)
            blank = 3;
            blankSize = mod(length(data),3);
        else
            blankSize = 0;
        end
        if blankSize ~= 0
            blankSize = blank-blankSize;
            blankData = zeros(blankSize,1);
            data1 = cat(1,data,blankData);
        else
            data1 = data;
        end
        [code, trellis, puncpat] = BCCEncode(data1,rate);

    %受信時の復号処理
    case 2
        trellis = cfgHARQ.trellis;
        puncpat = cfgHARQ.puncpat;
        code = BCCDecode(data,trellis,rate,puncpat);
        code = code(1:psdu,1);

    %再送信時の符号化処理
    case 3
        [code,~,puncpat] = BCCEncode(data,rate);
        code = code(psdu+1:psdu/rate,1);
        
    otherwise
        error('Unexpected HARQ mode. Only 1, 2, and 3 can use.')
    
end
end

