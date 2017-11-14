function [pscal] = pscal_lbd(PC, lbd_i, lbd_j)

matPC_sqnorm=diag(PC.PsiSqNorm);
pscal=(lbd_i)'*matPC_sqnorm*lbd_j;