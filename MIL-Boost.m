addpath(genpath('MIL'))
addpath(genpath('prtools'))
addpath(genpath('minfunc'))
%%
clear; clc; close all
prwarning(0)

load('XY_MIL.mat') %load data

rng(1)

fold_out=10;
fold_in=5;

global lrsel wlsel lropt wlopt;
lrsel=0;wlsel=0;lropt=0;wlopt=0;

lr=[0.0001 0.001 0.01 0.1 1 10 100];
wl=[5 10 15];
indice_opt_lr=[];
indice_opt_wl=[];

label_tot=[ones(size(Bag_past_negative,1),1); 2*ones(size(Bag_past_positive,1),1)];
Bag_tot=[Bag_past_negative; Bag_past_positive];

%%

%delete triglycerides
for i=1:size(Bag_tot,1)
   Bag_tot{i,1}(:,83)=[];
end
name_of_predictors(83)=[];

%delete glycaemia
for i=1:size(Bag_tot,1)
    Bag_tot{i,1}(:,82)=[];
end
name_of_predictors(82)=[];

%delete Hb1Ac
for i=1:size(Bag_tot,1)
   Bag_tot{i,1}(:,29)=[];
end
name_of_predictors(29)=[];


Bag_tot_aux=Bag_tot;
Bag_tot_aux=cell2mat(Bag_tot_aux);
a=sum(isnan(Bag_tot_aux));
th_nan=10; 
threshold=floor((th_nan*size(Bag_tot_aux,1))/100); 
idx=find((size(Bag_tot_aux,1)-a)<threshold);

Bag_tot_aux(:,idx)=[]; %delete features that contain missing values more than 90%
name_of_predictors(idx)=[];

%% Data imputation
Bag_tot_aux=knnimpute(Bag_tot_aux',1);
Bag_tot_aux=Bag_tot_aux.';
%%

seqq=[];
for u1=1:size(Bag_tot,1)
    seqq=[seqq; size(Bag_tot{u1,1},1)];
end

Bag_tot=mat2cell(Bag_tot_aux, seqq, size(Bag_tot_aux,2));

%%

Index_out=crossvalind('Kfold',label_tot,fold_out);

for out=1:fold_out
    disp(out)
    
    lab_training_out=label_tot(Index_out~=out);
    lab_test_out=label_tot(Index_out==out);
    
    selCVtrain=find(Index_out~=out);
    selCVtest=find(Index_out==out);
    training_out=[];
    for iu=1:numel(selCVtrain)
        training_out{iu,1}=Bag_tot{selCVtrain(iu),1};
    end
    test_out=[];
    for iu=1:numel(selCVtest)
        test_out{iu,1}=Bag_tot{selCVtest(iu),1};
    end
    
    PP=find(lab_training_out==2);
    NN=find(lab_training_out==1);
    
    BagTrainPos=[];
    BagTrainNeg=[];
    for iu=1:numel(PP)
        BagTrainPos{iu,1}=training_out{PP(iu),1};
    end
    
    for iu=1:numel(NN)
        BagTrainNeg{iu,1}=training_out{NN(iu),1};
    end
    BagTrainTot=[BagTrainNeg;BagTrainPos];
    traintot=cell2mat(BagTrainTot);
    [TRtot,mu_out,sigma_out]=zscore(traintot);
    
    seqp=[];
    seqn=[];
    for u1=1:size(BagTrainNeg,1)
        seqn=[seqn; size(BagTrainNeg{u1,1},1)];
    end
    for u2=1:size(BagTrainPos,1)
        seqp=[seqp; size(BagTrainPos{u2,1},1)];
    end
    sizecol=size(TRtot,2);
    sizerowP=size(cell2mat(BagTrainPos),1);
    sizerowN=size(cell2mat(BagTrainNeg),1);
    BagMatNeg=TRtot(1:sizerowN,:);
    BagMatPos=TRtot(sizerowN+1:end,:);
    BagTrainNegNorm=mat2cell(BagMatNeg, seqn, sizecol);
    BagTrainPosNorm=mat2cell(BagMatPos, seqp, sizecol);
    
    lab = genlab([sizerowN sizerowP]);
    idd=[];
    for iu=1:numel(BagTrainTot)
        idd(iu)=size(BagTrainTot{iu,1},1);
    end
   
    bagid=genlab(idd);
    milDatasetTR=genmil(TRtot,lab,bagid);
    milDataset2TR = positive_class(milDatasetTR,2);

    %% VALIDATION
    Index_in=crossvalind('Kfold',lab_training_out,fold_in);
    for in=1:fold_in
        disp(in)
        
        lab_training_in=lab_training_out(Index_in~=in);
        lab_test_in=lab_training_out(Index_in==in);
    
        selCVtrain_in=find(Index_in~=in);
        selCVtest_in=find(Index_in==in);
        training_in=[];
        for iu=1:numel(selCVtrain_in)
        training_in{iu,1}=BagTrainTot{selCVtrain_in(iu),1};
        end
        test_in=[];
        for iu=1:numel(selCVtest_in)
        test_in{iu,1}=BagTrainTot{selCVtest_in(iu),1};
        end
    
        PP_in=find(lab_training_in==2);
        NN_in=find(lab_training_in==1);
    
        BagTrainPos_in=[];
        BagTrainNeg_in=[];
        for iu=1:numel(PP_in)
        BagTrainPos_in{iu,1}=training_in{PP_in(iu),1};
        end
    
        for iu=1:numel(NN_in)
        BagTrainNeg_in{iu,1}=training_in{NN_in(iu),1};
        end
        BagTrainTot_in=[BagTrainNeg_in;BagTrainPos_in];
        traintot_in=cell2mat(BagTrainTot_in);
        [TRtot_in,mu_in,sigma_in]=zscore(traintot_in);
    
        seqp_in=[];
        seqn_in=[];
        for u1=1:size(BagTrainNeg_in,1)
        seqn_in=[seqn_in; size(BagTrainNeg_in{u1,1},1)];
        end
        for u2=1:size(BagTrainPos_in,1)
        seqp_in=[seqp_in; size(BagTrainPos_in{u2,1},1)];
        end
        sizecol_in=size(TRtot_in,2);
        sizerowP_in=size(cell2mat(BagTrainPos_in),1);
        sizerowN_in=size(cell2mat(BagTrainNeg_in),1);
        BagMatNeg_in=TRtot_in(1:sizerowN_in,:);
        BagMatPos_in=TRtot_in(sizerowN_in+1:end,:);
        BagTrainNegNorm_in=mat2cell(BagMatNeg_in, seqn_in, sizecol_in);
        BagTrainPosNorm_in=mat2cell(BagMatPos_in, seqp_in, sizecol_in);
    
        lab_in = genlab([sizerowN_in sizerowP_in]);
        idd_in=[];
        for iu=1:numel(BagTrainTot_in)
        idd_in(iu)=size(BagTrainTot_in{iu,1},1);
        end
   
        bagid_in=genlab(idd_in);
        milDatasetTR_in=genmil(TRtot_in,lab_in,bagid_in);
        milDataset2TR_in = positive_class(milDatasetTR_in,2);
        
        testtot_in=cell2mat(test_in);
        C_in = bsxfun(@minus, testtot_in, mu_in);
        sigma_in(sigma_in==0)=eps;
        TEtot_in = bsxfun(@rdivide, C_in, sigma_in);
        seqtest_in=[];
        for u1=1:size(test_in,1)
        seqtest_in=[seqtest_in; size(test_in{u1,1},1)];
        end
     
        BagTEnorm_in=mat2cell(TEtot_in,seqtest_in,sizecol_in);
        idpos_in=find(lab_test_in==2);
        idneg_in=find(lab_test_in==1);
        BagTENegNorm_in=[];
        BagTEPosNorm_in=[];
        for jj=1:numel(idneg_in)
        BagTENegNorm_in{jj,1}=BagTEnorm_in{idneg_in(jj),1};
        end
        for jj=1:numel(idpos_in)
        BagTEPosNorm_in{jj,1}=BagTEnorm_in{idpos_in(jj),1};
        end
        sizerowPt_in=size(cell2mat(BagTEPosNorm_in),1);
        sizerowNt_in=size(cell2mat(BagTENegNorm_in),1);
        TEtot_in=[cell2mat(BagTENegNorm_in); cell2mat(BagTEPosNorm_in)];
   
        labt_in = genlab([sizerowNt_in sizerowPt_in]);
        iddT_in=[];
        for iu=1:numel(BagTEnorm_in)
        iddT_in(iu)=size(BagTEnorm_in{iu,1},1);
        end
        bagid_in=genlab(iddT_in);
    
        milDatasetTE_in=genmil(TEtot_in,labt_in,bagid_in);
        milDataset2TE_in = positive_class(milDatasetTE_in,2);
    
        for h=1:length(lr)
            for hh=1:length(wl)
            setGlobalx(lr(h),wl(hh),0,0);
            w_opt=boosting_mil_in(milDataset2TR_in);
    
            yp_opt2=[];
            yp_opt=milDataset2TE_in*w_opt*labeld;
    
            for jj=1:numel(lab_test_in)
                if yp_opt(jj)=='p'
                    yp_opt2(jj,1)=2;
                else
                    yp_opt2(jj,1)=1;
                end
            end
    
            accuacy_opt{in}(h,hh)=(sum(yp_opt2==lab_test_in))/length(lab_test_in);
           [macro_opt{in}(h,hh), precision_opt{in}(h,hh), recall_opt{in}(h,hh)] = my_micro_macro(yp_opt2,lab_test_in);
            posterior_in=double(milDataset2TE_in*w_opt); %positive and negative  
           [~,~,~,AUC_opt{in}(h,hh)] = perfcurve(lab_test_in,posterior_in(:,1),2); 
    
            end
        end
    end
    
    recall_opt_mean=zeros(length(lr), length(wl));
    
    for tt=1:in
    recall_opt_mean=recall_opt_mean+recall_opt{tt};
    end
    
    recall_opt_mean=recall_opt_mean/fold_in;
    
    [v,l]=max(recall_opt_mean(:));
    
    [optlr, optwl]=ind2sub(size(recall_opt_mean),l);
    
    indice_opt_lr(out)=optlr;
    indice_opt_wl(out)=optwl;
    
    %% TRAIN
    
    testtot=cell2mat(test_out);
    C_out = bsxfun(@minus, testtot, mu_out);
    sigma_out(sigma_out==0)=eps;
    TEtot = bsxfun(@rdivide, C_out, sigma_out);
    seqtest=[];
    for u1=1:size(test_out,1)
        seqtest=[seqtest; size(test_out{u1,1},1)];
    end
     
    BagTEnorm=mat2cell(TEtot,seqtest,sizecol);
    idpos=find(lab_test_out==2);
    idneg=find(lab_test_out==1);
    BagTENegNorm=[];
    BagTEPosNorm=[];
    for jj=1:numel(idneg)
        BagTENegNorm{jj,1}=BagTEnorm{idneg(jj),1};
    end
    for jj=1:numel(idpos)
        BagTEPosNorm{jj,1}=BagTEnorm{idpos(jj),1};
    end
    sizerowPt=size(cell2mat(BagTEPosNorm),1);
    sizerowNt=size(cell2mat(BagTENegNorm),1);
    TEtot=[cell2mat(BagTENegNorm); cell2mat(BagTEPosNorm)];
   
    labt = genlab([sizerowNt sizerowPt]);
    iddT=[];
    for iu=1:numel(BagTEnorm)
        iddT(iu)=size(BagTEnorm{iu,1},1);
    end
    bagid=genlab(iddT);
    
    milDatasetTE=genmil(TEtot,labt,bagid);
    milDataset2TE = positive_class(milDatasetTE,2);
    
    setGlobalx(0,0,lr(optlr),wl(optwl));
    w=boosting_mil_out(milDataset2TR); % Number of weak learners 50

    %% TEST   
    
    yp2=[];
    yp=milDataset2TE*w*labeld;
    
    for jj=1:numel(lab_test_out)
        if yp(jj)=='p'
            yp2(jj,1)=2;
        else
            yp2(jj,1)=1;
        end
    end
    
    accDDtest(out,1)=(sum(yp2==lab_test_out))/length(lab_test_out);
    ConfDDtest{out,1}=confusionmat(lab_test_out,yp2);
    [FDDmacro(out,1), FDDprecision(out,1), FDDrecall(out,1)] = my_micro_macro(yp2,lab_test_out);
    
    posterior=double(milDataset2TE*w); %positive and negative  
    [X_ext{out},Y_ext{out},T_ext{out},FDDAUC(out)] = perfcurve(lab_test_out,posterior(:,1),2);
end

ConfDDtesttot=zeros(2,2);

for out=1:fold_out
    ConfDDtesttot=ConfDDtesttot+ConfDDtest{out,1};
end

ConfDDtesttot=(ConfDDtesttot./repmat(sum(ConfDDtesttot,2),1,2)*100);

FDDaccuracy_tot=mean(accDDtest);
FDDmacro_tot=mean(FDDmacro);
FDDprecision_tot=mean(FDDprecision);
FDDrecall_tot=mean(FDDrecall);
FDDAUC_tot=mean(FDDAUC);

save('resultsMILBoost_knn_noTyg_10')
