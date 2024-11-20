function err = vouerror(rep,y)

err = rep > 0;   % DARE failed: show-stopper!

if nargin < 2 || isempty(y)
    if err
        fprintf(2,'CARE ERROR: ');
        switch rep
            case 1, fprintf(2,'solution accuracy is poor\n');
            case 2, fprintf(2,'solution is not finite\n');
            case 3, fprintf(2,'eigenvalues on the imaginary axis\n');
            case 4, fprintf(2,'pencil is singular\n');
        end
        return
    end
else
    if err
        fprintf(2,'CARE ERROR for source %d: ',y);
        switch rep
            case 1, fprintf(2,'solution accuracy is poor\n');
            case 2, fprintf(2,'solution is not finite\n');
            case 3, fprintf(2,'eigenvalues on the imaginary axis\n');
            case 4, fprintf(2,'pencil is singular\n');
        end
        return
    end
end
