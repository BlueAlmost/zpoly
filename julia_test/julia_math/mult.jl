using DSP

a = [1.1, 2.2, -3.3];
b = [1.1, 2.2, 3.3, -4.4];

c = conv(a,b);

aa = [1.1, 2.2, -3.3] + im*[1.0, 3.0, 4.0];
bb = [1.1, 2.2, 3.3, -4.4] + im*[2.0, 1.0, 3.0, -3.0];

cc = conv(aa,bb);


