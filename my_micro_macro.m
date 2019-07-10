function [ fmacro, macroprecision, macrorecall ] = my_micro_macro( pred_label, orig_label)
%computer micro and macro: precision, recall and fscore
%Sandy wltongxing@163.com
%micro>macro?
mat=confusionmat(orig_label, pred_label);
%label_unique=unique([orig_label(:);pred_label(:)]);
%     microTP=0;
%     microFP=0;
%     microFN=0;
len=size(mat,1);
macroTP=zeros(len,1);
macroFP=zeros(len,1);
macroFN=zeros(len,1);
macroP=zeros(len,1);
macroR=zeros(len,1);
macroF=zeros(len,1);
for i=1:len
    macroTP(i)=mat(i,i);
    macroFP(i)=sum(mat(:, i))-mat(i,i);
    macroFN(i)=sum(mat(i,:))-mat(i,i);
    if macroTP(i)==0 && macroFP(i)==0 
        macroP(i)=0;
    else    
    macroP(i)=macroTP(i)/(macroTP(i)+macroFP(i)+eps);
    end
    if macroTP(i)==0 && macroFN(i)==0 
        macroR(i)=0;
    else    
    macroR(i)=macroTP(i)/(macroTP(i)+macroFN(i)+eps);
    end
    if macroP(i)==0 && macroR(i)==0 
        macroF(i)=0;
    else    
    macroF(i)=2*macroP(i)*macroR(i)/(macroP(i)+macroR(i));
    end
end
% macro.precision=mean(macroP);
% macro.recall=mean(macroR);
fmacro=mean(macroF);
macroprecision=mean(macroP);
macrorecall=mean(macroR);

micro.precision=sum(macroTP)/(sum(macroTP)+sum(macroFP));
micro.recall=sum(macroTP)/(sum(macroTP)+sum(macroFN));
micro.fscore=2*micro.precision*micro.recall/(micro.precision+micro.recall);
end