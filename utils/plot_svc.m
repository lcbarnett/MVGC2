function plot_svc(sval,svc,mosvc,rmax,plotm)

% Plot Bauer's Singular Value Criterion (SVC) for state-space subspace
% SS model estimation.

mo = (1:rmax)';

if mosvc == rmax, wsvc = '*'; else wsvc = ''; end

gap = 0.05;
ssvc = gap+(1-gap)*(svc-min(svc))/(max(svc)-min(svc));

if ischar(plotm) % Gnuplot

	gpname = 'sssvc';
	gpstem = fullfile(tempdir,gpname);
	gpdat = [mo ssvc sval];
	gp_write(gpstem,gpdat);

	gp = gp_open(gpstem,plotm,[Inf,0.6]);

	fprintf(gp,'datfile = "%s.dat"\n',gpname);

	fprintf(gp,'\nset grid lt 1 lc "dark-grey"\n');
	fprintf(gp,'set xr[0:%g]\n',rmax);
	fprintf(gp,'set xlabel "State space dimension"\n');
	fprintf(gp,'set arrow 1 from first %g,graph 0 to first %g,graph 1 lt 1 lc "red" lw 2 nohead\n',mosvc,mosvc);

	fprintf(gp,'\nset multiplot title "SS SVC model order selection (CCA, max = %d)\\\n" layout 2,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n',rmax);

	fprintf(gp,'\nset title "Singular Value Criterion (SVC)"\n');
	fprintf(gp,'set ytics 0.2\n');
	fprintf(gp,'set yr[0:1.05]\n');
	fprintf(gp,'set ylabel "SVC (scaled)"\n');
	fprintf(gp,'set key top right Left rev\n');
	fprintf(gp,'plot \\\n');
	fprintf(gp,'datfile u 1:2 w linespoints pt 6 ps 1.4 t "SVC (opt = %2d%c)"\n',mosvc,wsvc);

	fprintf(gp,'\nset title "Singular values"\n');
	fprintf(gp,'unset logs y\n');
	fprintf(gp,'set ytics auto format ''%% h''\n');
	fprintf(gp,'set yr [0:*]\n');
	fprintf(gp,'set ytics 0.2 nomirror\n');
	fprintf(gp,'set ylabel "Singular value"\n');
	fprintf(gp,'plot datfile u 1:3 w boxes fs solid 0.25 not\n');

	fprintf(gp,'\nunset multiplot\n');

	gp_close(gp,gpstem,plotm);

else % Matlab

	if plotm == 0, figure; else, figure(plotm); end; clf;

	xlims = [0 rmax];

	subplot(2,1,1);
	plot(mo,ssvc,'o-');
	grid on
	title('Singular Value Criterion (SVC)');
	ylabel('SVC (scaled)');
	xlabel('State space dimension');
	legend(sprintf('SVC (opt = %d%c)',mosvc,wsvc));
	xlim(xlims);
	ylim([0 1+gap]);
	xline(mosvc,'r','HandleVisibility','off');

	subplot(2,1,2); % SVC
	bar(mo,sval,1.01,'FaceColor',[0.65 0.75 1]);
	grid on
	xlim(xlims);
	xlabel('State space dimension');
	ylabel('Singular value');
	title('Singular values');
	xline(mosvc,'r','HandleVisibility','off');

	axes('Units','Normal');
	h = title(sprintf('SS SVC model order selection (CCA, max = %d)\n\n',rmax),'FontSize',13);
	set(gca,'visible','off')
	set(h,'visible','on')
end
