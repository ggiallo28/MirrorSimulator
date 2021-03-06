function [fitresult, gof] = createFit(yy, aa)
%CREATEFIT(YY,AA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : yy
%      Y Output: aa
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 01-Dec-2016 09:32:26


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( yy, aa );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );


