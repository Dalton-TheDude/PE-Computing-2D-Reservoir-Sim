syms p
y=p*(p-5000);
b=1.2/(1+y+(0.5*(y^2)));
f(p)=diff(b,p);
f(5)
