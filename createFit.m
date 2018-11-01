function [fitresult, gof] = createFit(u, lgcv)
%CREATEFIT(U,LGCV)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : u
%      Y Output: lgcv
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 01-Nov-2018 00:21:37


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( u, lgcv );

% Set up fittype and options.
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.994226483476193 0.843260710269445];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'lgcv vs. u', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel u
ylabel lgcv
grid on


