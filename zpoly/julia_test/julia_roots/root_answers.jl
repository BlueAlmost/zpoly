using Polynomials

ar = [1.0, 1.2, -3.2, 2.1, 3.3, 3.4, -4.1];
PR = Polynomial(ar);
Rts_R = roots(PR)

println(" ")
println(" ")
println("real polynomial")
for i = 1:length(ar)-1
   println("   try std.testing.expectApproxEqAbs(x[",i-1,"].re,  ", real(Rts_R[i]), ", eps);")
end
println(" ")
for i = 1:length(ar)-1
   println("   try std.testing.expectApproxEqAbs(x[",i-1,"].im,  ", imag(Rts_R[i]), ", eps);")
end

ac = [1.0, 1.2, -3.2, 2.1, 3.3, 3.4, -4.1] + im*[-1.1, 1.12, 2.12, 1.0, 0.13, -1.1, -1.2] 
PC = Polynomial(ac);
Rts_C = roots(PC)

println(" ")
println(" ")
println("complex polynomial")
for i = 1:length(ar)-1
   println("   try std.testing.expectApproxEqAbs(x[",i-1,"].re,  ", real(Rts_C[i]), ", eps);")
end
println(" ")
for i = 1:length(ar)-1
   println("   try std.testing.expectApproxEqAbs(x[",i-1,"].im,  ", imag(Rts_C[i]), ", eps);")
end

