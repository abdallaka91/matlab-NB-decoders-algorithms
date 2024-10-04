function [HD1, HD1_complx] = HD_QAM(y, constl1)
% y1 = bsxfun(@minus,y, constl1.');
y1 = y- constl1.';

[~,HD1] = min(abs(y1),[],2);
HD1_complx = constl1(HD1);
HD1 = HD1-1;