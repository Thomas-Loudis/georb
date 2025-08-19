
mjd0 = 53340

mjd(1,1) = mjd0 - 1;
mjd(2,1) = mjd0;
mjd(3,1) = mjd0 + 1;
mjd(4,1) = mjd0 + 2;

[eop] = eop_read(filename,mjd)

X_matrix = eop(:,1)
Y_matrix = eop(:,2)
xint  = mjd0 + 100 / 86400
dpint = 4
[yint] = interp_Lagrange(X_matrix,Y_matrix,xint,dpint)
