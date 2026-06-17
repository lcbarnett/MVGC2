function [cvtab,sigvals] = jotest_cvtab(model)

% Johansen test critical value tables and significance levels

CV = load('Data_JCITest');
switch model
	case 'H2',  cvtab = squeeze(CV.JCV(:,:,1,:));
	case 'H1*', cvtab = squeeze(CV.JCV(:,:,2,:));
	case 'H1',  cvtab = squeeze(CV.JCV(:,:,3,:));
	case 'H*',  cvtab = squeeze(CV.JCV(:,:,4,:));
	case 'H',   cvtab = squeeze(CV.JCV(:,:,5,:));
	otherwise,  error('Unknown cointegration model');
end

sigvals = [0.001 0.005:0.005:0.10 0.125:0.025:0.20 0.80:0.025:0.875 0.90:0.005:0.995 0.999];
