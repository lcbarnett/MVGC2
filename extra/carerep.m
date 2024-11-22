function err = carerep(errcode)

if errcode ~= 0
	switch errcode
		case 1, err = 'CARE solution accuracy is poor';
		case 2, err = 'CARE solution is not finite';
		case 3, err = 'CARE eigenvalues on the imaginary axis';
		case 4, err = 'CARE pencil is singular';
		otherwise, error('unknown CARE error');
	end
else
	err = '';
end
