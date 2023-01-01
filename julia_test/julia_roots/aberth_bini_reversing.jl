using Polynomials
using PyPlot

function upper_convex_hull(k, a)
   la = log2.(abs.(a));
   L_k = [ k[1], k[2] ]
   L_la = [ la[1], la[2] ]
   for i = 3:length(a)
      append!(L_k, k[i])
      append!(L_la, la[i])
      while ( length(L_la)>2 )
         ab_x = L_k[end-1]-L_k[end-2]
         ab_y = L_la[end-1]-L_la[end-2]
         bc_x = L_k[end]-L_k[end-1]
         bc_y = L_la[end]-L_la[end-1]
         cross = ab_x*bc_y - ab_y*bc_x
         if( cross >= 0)
            deleteat!(L_k, length(L_k)-1)
            deleteat!(L_la, length(L_la)-1)
         else
            break
         end
      end
   end
   return L_k, L_la
end

function compute_modulii(k, a)
   q = length(k)
   u = zeros(q-1)
   for i = 2:q
      u[i-1] = abs( a[k[i-1]+1] / a[k[i]+1]) ^ (1 / (k[i] - k[i-1]))
   end
   return u
end

function initial_vals(k, u)
   q = length(k)
   n = k[end]
   x = zeros(k[end]) + im*zeros(k[end])
   sigma = 0.7;

   println(" ")
   for i = 1:q-1
      for j = 1: k[i+1]-k[i]
         alpha = (2*pi*j)/(k[i+1]-k[i])
         beta = (2*pi*i/n)
         phi = exp(im* (alpha + beta + sigma))
         x[ k[i]+j] = u[ i ] * phi
      end
   end
   println(" ")
   return x
end

#-----------------------------------------------------------------

a = [1.0, 1.2, -3.2, 2.1, 3.3, 3.4, -4.1];
# a = [1.0, 6.5, -2.3433, 7.34534, 1.2143, -3.11112, 2.1, 6.11, 3.3, 3.4, -4.1];
n = length(a)-1;

global sx = 0
mu = 2^-53
for i=0:n
   global sx
   sx = sx + abs(a[i+1])*(1+4*i)
end

P = Polynomial(a);
R = roots(P)

# set initial approximation values -----------------
fig = figure()
fig.set_figwidth(8)
fig.set_figheight(8)

# compute convex hull, keep absissa values (k) - there are "q" of them
stem( collect(0:n), log2.(abs.(a)))
k, hull = upper_convex_hull( collect(0:n), a)

# figure()
plot(k, hull, color = "orange")
pause(1)

figure()
clf()

# compute modulii u
u = compute_modulii(k, a)
println("modulii: ", u)
unit_circle = exp.(im*2*pi*collect(0:200)/200)
plot( real(unit_circle), imag(unit_circle))

for i=1:length(u)
   plot( u[i]*real(unit_circle), u[i]*imag(unit_circle), ls=":")
end


# initializations
x = initial_vals(k, u)

xinits = x

plot( real(x), imag(x), marker="x", ls="none")
pause(1)


for iter = 1:50
   for i = 1:n

      if(abs(x[i]) <= 1)

         px = evalpoly(x[i], a)
         d = collect(1:n).*a[2:end];
         dx = evalpoly(x[i], d)

         pxdx= px/dx
      else
         # println("reversing")
         b = reverse(a)
         y = 1/x[i]
         py = evalpoly(y, b)
         d = collect(1:n).*b[2:end];
         dy = evalpoly(y, d)

         pxdx = 1/(n*y - y*y*(dy/py))

      end

      sum = 0;
      for j = 1 : n
         if (j != i)
            sum = sum + (1/(x[i]-x[j]))
         end
      end
      x[i] = x[i] - (pxdx / (1-pxdx*sum))

   end

   plot( real(x), imag(x), marker=".", ls="none")
   pause(0.1)

good_enough = false
   for i = 1:n
      if( abs(evalpoly(x[i],a)) < mu*sx )
         good_enough = true
         # println("p(x[", i, "] = ", abs(evalpoly(x[i],a)), "   good_enough")
      else
         good_enough = false
      end
   end
   # println(" ")

   if( good_enough)
      # println(iter, " iterations performed")
      # break
   end
end


global mine = 0
global jul = 0
for i = 1:n
   global mine = mine + abs(evalpoly(x[i], a))
   global jul = jul + abs(evalpoly(R[i], a))
end
println("mine: ", mine)
println(" jul: ", jul)




println(" ")
println("polynomial roots")
for i = 1:length(x)
   println("   try std.testing.expectApproxEqAbs(x[",i-1,"].re,  ", real(x[i]), ", eps);")
end
println(" ")
for i = 1:length(x)
   println("   try std.testing.expectApproxEqAbs(x[",i-1,"].im,  ", imag(x[i]), ", eps);")
end



