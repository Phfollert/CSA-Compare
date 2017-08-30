csa = "sada";
op = "extract";
folder = strcat("/home/nacho/Desktop/sdsl-lite/CSA COMPARE/",csa,"/");
a = csvread(strcat(folder,"dna50mb_",op,".txt"));
b = csvread(strcat(folder,"dna100mb_",op,".txt"));
c = csvread(strcat(folder,"dna200mb_",op,".txt"));
d = csvread(strcat(folder,"english50mb_",op,".txt"));
e = csvread(strcat(folder,"english100mb_",op,".txt"));
f = csvread(strcat(folder,"english200mb_",op,".txt"));
ae = csvread(strcat(folder,"dna50mb_",op,"_error.txt"));
be = csvread(strcat(folder,"dna100mb_",op,"_error.txt"));
ce = csvread(strcat(folder,"dna200mb_",op,"_error.txt"));
de = csvread(strcat(folder,"english50mb_",op,"_error.txt"));
ee = csvread(strcat(folder,"english100mb_",op,"_error.txt"));
fe = csvread(strcat(folder,"english200mb_",op,"_error.txt"));

x = 2000:1:2099;
plot(x,a,".-b",x,b,".-g",x,c,".-y",x,d,".-m",x,e,".-r", x,f,".-k")
title (cstrcat(csa,": Gráfico tiempo vs largo de prefijo para operación ",op));
xlabel ("Largo del prefijo");
ylabel ("Tiempo en nanosegundos");
legend ("50mb dna", "100mb dna", "200mb dna", "50mb english", "100mb english", "200mb english");