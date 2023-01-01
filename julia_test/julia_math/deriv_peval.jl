# DERIV EVAL --------------------------------------------
xr = -4.4
xc = -0.2 + im*0.1;

ar = [ 1.1, 2.2, -3.3, 1.2];
ac = [ 1.2, -0.11, 1.20, -0.40] + im*[2.31, 2.70, -1.33, 2.34];


# FORWARD
# real, real
deriv_rr = ar[2] + 2*ar[3]*xr + 3*ar[4]*xr^2

# real, complex
deriv_rc = ar[2] + 2*ar[3]*xc + 3*ar[4]*xc^2

# real, complex
deriv_cr = ac[2] + 2*ac[3]*xr + 3*ac[4]*xr^2

# complex, complex
deriv_cc = ac[2] + 2*ac[3]*xc + 3*ac[4]*xc^2


# REVERSE
br = reverse(ar);
bc = reverse(ac);
# real, real
deriv_rev_rr = br[2] + 2*br[3]*xr + 3*br[4]*xr^2

# real, complex
deriv_rev_rc = br[2] + 2*br[3]*xc + 3*br[4]*xc^2

# real, complex
deriv_rev_cr = bc[2] + 2*bc[3]*xr + 3*bc[4]*xr^2

# complex, complex
deriv_rev_cc = bc[2] + 2*bc[3]*xc + 3*bc[4]*xc^2

