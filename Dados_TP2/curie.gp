# Define the data file
datafile = 'X&invX_T__GD2O3.txt'

# Define the Curie-Weiss fitting function
#f(x) = 1/(C/(x - Tp) + D)
f(x) = (C/(x - Tp)) + D

# Provide initial guesses for the parameters
C = 1e-3
Tp = 1e-3
D = 1e-3

# Perform the fit
fit [100:300] f(x) datafile u 1:2 via C, Tp, D

# Print the values of the parameters
print sprintf("C = %g +/- %g, Tp = %g +/- %g, D = %g +/- %g", C, C_err, Tp, Tp_err, D, D_err)

miu = (C*4)**0.5
delta_miu = 0.5*0.5*((C*0.5)**(-0.5)) * C_err

print sprintf("miu = %g +/- %g", miu, delta_miu)

# Plot the data and the fitted function
set title "1/X_{molar} vs Temperature"
set xlabel "Temperature (K)"
set ylabel "X_{molar}"
plot datafile u 1:2 with points title 'Data', f(x) with lines title 'Fit'

# Wait for user to close the interactive plot
pause -1 "Hit any key to continue"

# Save the plot to a file
set terminal png
set output 'fit_plot.png'
plot datafile u 1:2 with points title 'Data', f(x) with lines title 'Fit'
set output
