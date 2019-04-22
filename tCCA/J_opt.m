function [x, fval] = J_opt(x, CORR, MSE, PVAL, FSCORE, CHNO, fact)
% J_OPT Cost function to optimize
% uses results from run_CCA 
% INPUTS
% CORR, MSE, PVAL, FSCORE, CHNO are 4D matrices with results from run_CCA 
% averaged across subjects, folds and channels: 
%   (HbO/R (2) x tlags (11) x stpsize (12) x cthresh (10)) 

% x: input vector with indices [tlag, stsize, cthresh]
% fact: struct with factors (weights) for adapting J
%   .corr
%   .mse
%   .pval
%   .fscore
%   .chno
%   .HbO
%   .HbR

% Input data dimensions: 
% #SBJ x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres

fval = fact.HbO*fact.corr*CORR(1,x(1),x(2),x(3)) ...
    + fact.HbR*fact.corr*CORR(2,x(1),x(2),x(3)) ...
    + fact.HbO*fact.mse*MSE(2,x(1),x(2),x(3)) ...
    + fact.HbR*fact.mse*MSE(2,x(1),x(2),x(3)) ...
    + fact.HbO*fact.pval*PVAL(2,x(1),x(2),x(3)) ...
    + fact.HbR*fact.pval*PVAL(2,x(1),x(2),x(3)) ...





end

