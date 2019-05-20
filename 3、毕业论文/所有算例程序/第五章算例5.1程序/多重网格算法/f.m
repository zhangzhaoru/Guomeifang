function f1=f(x,y,t,alpha1,alpha2,Kx,Ky,Ca1,Ca2)
gx1=(24/gamma(5-alpha1))*x.^(4-alpha1)-(12/gamma(4-alpha1))*x.^(3-alpha1)+(2/gamma(3-alpha1))*x.^(2-alpha1);
gx2=(24/gamma(5-alpha1))*(1-x).^(4-alpha1)-(12/gamma(4-alpha1))*(1-x).^(3-alpha1)+(2/gamma(3-alpha1))*(1-x).^(2-alpha1);
gy1=(24/gamma(5-alpha2))*y.^(4-alpha2)-(12/gamma(4-alpha2))*y.^(3-alpha2)+(2/gamma(3-alpha2))*y.^(2-alpha2);
gy2=(24/gamma(5-alpha2))*(1-y).^(4-alpha2)-(12/gamma(4-alpha2))*(1-y).^(3-alpha2)+(2/gamma(3-alpha2))*(1-y).^(2-alpha2);
f1=-10*exp(-t)*x.^2.*(1-x).^2.*y.^2.*(1-y).^2+...
    100*exp(-2*t)*x.^4.*(1-x).^4.*y.^4.*(1-y).^4+...
    10*Kx*Ca1*exp(-t)*y.^2.*(1-y).^2.*(-gx1-gx2)+...
    10*Ky*Ca2*exp(-t)*x.^2.*(1-x).^2.*(-gy1-gy2);
end

