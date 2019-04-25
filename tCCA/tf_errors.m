% True Positive *Rate* and False Positive *Rate*
TPR_SS = Ch_TP_SS./(Ch_TP_SS + Ch_FN_SS);
FPR_SS = Ch_FP_SS./(Ch_FP_SS + Ch_TN_SS);

fPR_CCA = Ch_TP_CCA./(Ch_TP_CCA + Ch_FN_CCA);
FPR_CCA = Ch_FP_CCA./(Ch_FP_CCA + Ch_TN_CCA);

% get F-score
% SS & CCA
Precision_SS = Ch_TP_SS ./(Ch_TP_SS + Ch_FP_SS);
Recall_SS = Ch_TP_SS ./(Ch_TP_SS + Ch_FN_SS);
F_score_SS = 2 * (Precision_SS .* Recall_SS)./(Precision_SS + Recall_SS);

Precision_CCA = Ch_TP_CCA ./(Ch_TP_CCA + Ch_FP_CCA);
Recall_CCA = Ch_TP_CCA ./(Ch_TP_CCA + Ch_FN_CCA);
F_score_CCA = 2 * (Precision_CCA .* Recall_CCA)./(Precision_CCA + Recall_CCA);

% reshape all to "# of sbjs x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
Ch_TP_SS = reshape(Ch_TP_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
Ch_FP_SS = reshape(Ch_FP_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
Ch_FN_SS = reshape(Ch_FN_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
Ch_TN_SS = reshape(Ch_TN_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
F_score_SS = reshape(F_score_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
TPR_SS = reshape(TPR_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
FPR_SS = reshape(FPR_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));

Ch_TP_CCA = reshape(Ch_TP_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
Ch_FP_CCA = reshape(Ch_FP_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
Ch_FN_CCA = reshape(Ch_FN_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
Ch_TN_CCA = reshape(Ch_TN_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
F_score_CCA = reshape(F_score_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
fPR_CCA = reshape(fPR_CCA, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
FPR_CCA = reshape(FPR_CCA, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
