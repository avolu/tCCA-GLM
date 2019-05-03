%% # of TP/FP/FN/TN channels
% for SS and CCA methods:
% dimensions of DET_CCA
% # Subjects x #CH x 2(Hbo+HbR) x 2 (cv split) x AUX combinations
foo_CCA = permute(DET_CCA,[2 1 3 4 5]); % swapping to get a matrix with ch X everything else
% #CH x # Subjects x 2(Hbo+HbR) x 2 (cv split) x AUX combinations
foo_CCA = reshape(foo_CCA, size(foo_CCA,1), size(foo_CCA,2)*size(foo_CCA,3)*size(foo_CCA,4)*size(foo_CCA,5));
% ROCLAB.name = {'TP','FP','FN','TN', 'PRND'};
for i = 1:size(foo_CCA,2)
    % CCA
    Ch_TP_CCA(i) = sum(foo_CCA(:,i)==1);
    Ch_FP_CCA(i) = sum(foo_CCA(:,i)==-1);
    Ch_FN_CCA(i) = sum(foo_CCA(:,i)==2);
    Ch_TN_CCA(i) = sum(foo_CCA(:,i)==-2);
end

% True Positive *Rate* and False Positive *Rate*
TPR_CCA = Ch_TP_CCA./(Ch_TP_CCA + Ch_FN_CCA);
FPR_CCA = Ch_FP_CCA./(Ch_FP_CCA + Ch_TN_CCA);

% get F-score
%  CCA

Precision_CCA = Ch_TP_CCA ./(Ch_TP_CCA + Ch_FP_CCA);
Recall_CCA = Ch_TP_CCA ./(Ch_TP_CCA + Ch_FN_CCA);
F_score_CCA = 2 * (Precision_CCA .* Recall_CCA)./(Precision_CCA + Recall_CCA);

% reshape all to "# of sbjs x 2(Hbo+HbR) x 2 (cv split) x AUX combinations
Ch_TP_CCA = reshape(Ch_TP_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
Ch_FP_CCA = reshape(Ch_FP_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
Ch_FN_CCA = reshape(Ch_FN_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
Ch_TN_CCA = reshape(Ch_TN_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
F_score_CCA = reshape(F_score_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
TPR_CCA = reshape(TPR_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
FPR_CCA = reshape(FPR_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5));
