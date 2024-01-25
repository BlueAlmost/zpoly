using Polynomials

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


P = Polynomial(a);
R = roots(P)

# set initial approximation values -----------------

# compute convex hull, keep absissa values (k) - there are "q" of them
k, hull = upper_convex_hull( collect(0:n), a)

# compute modulii u
u = compute_modulii(k, a)
println("modulii: ", u)


# initializations
x = initial_vals(k, u)

xinits = x

MAX_ITER = 11

println(" ")
for iter = 1:MAX_ITER
   for i = 1:n

      px = evalpoly(x[i], a)
      # println("     i: ", i, "   ", px)

      d = collect(1:n).*a[2:end];
      dx = evalpoly(x[i], d)
      # println("     i: ", i, "   ", dx)

      pxdx= px/dx
      # println("     i: ", i, "   ", pxdx)

      sum = 0;
      for j = 1 : n
         if (j != i)

            # println("             diff_x: ", (1/(x[i]-x[j])))
            sum = sum + (1/(x[i]-x[j]))
         end
      end
      
      # println("     i: ", i, "   ", sum)

      # println("  update: ", pxdx /(1-pxdx*sum))
      x[i] = x[i] - (pxdx / (1-pxdx*sum))

   end
end
println(" ")


global mine = 0
global jul = 0
println(" ")
for i = 1:n

   println("x[", i-1, "]: ( ", real(x[i]), " , ", imag(x[i]), " )")

   global mine = mine + abs(evalpoly(x[i], a))
   global jul = jul + abs(evalpoly(R[i], a))
end

println(" ")
println(" ")
println("mine: ", mine)
println(" jul: ", jul)
