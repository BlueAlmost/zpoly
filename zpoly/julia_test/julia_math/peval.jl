# REGULAR EVAL --------------------------------------------
# real, real
a = [ 1.1, 2.2, -3.3, 1.2];
x = -4.4
frr = evalpoly(x,a)

# real, complex
xc = -4.4 + im*1.0;
frc = evalpoly(xc,a)

# complex, real
ac = [ 1.2, -0.11, 1.20, -0.40] + im*[2.31, 2.70, -1.33, 2.34];
fcr = evalpoly(x,ac)

# complex, complex
ac = [ 1.2, -0.11, 1.20, -0.40] + im*[2.31, 2.70, -1.33, 2.34];
xc = -0.2 + im*0.1
fcc = evalpoly(xc,ac)


# REVERSE EVAL --------------------------------------
b = reverse(a)
rev_frr = evalpoly(x, b);

bc = reverse(ac)
rev_frc = evalpoly(xc, b);

bc = reverse(ac)
rev_fcr = evalpoly(x, bc);

bc = reverse(ac)
rev_fcc = evalpoly(xc, bc);


# df = a[2] + 2*a[3]*x + 3*a[4]*x^2
# db = b[2] + 2*b[3]*x + 3*b[4]*x^2
# rdf = a[3] + 2*a[2]*x + 3*a[1]*x^2
# ddff = aa[2] + 2*aa[3]*xx + 3*aa[4]*xx^2
# rrddff = aa[3] + 2*aa[2]*xx + 3*aa[1]*xx^2
