y = -10:10;
x = 0:150;
k1 = 2;
k2 = 1;
p_inf = 0.2;
xd = 70;
x_bar = x/xd;

c2 = (1 - p_inf)/(2^(k1+1)*(k2/k1) - 1 - tanh(k2));
c1 = 1 - p_inf + c2*(1 + tanh(k2));

px = c1*(1+x_bar).^(-k1) + c2*tanh(k2*(x_bar-1)) - c2 + p_inf;
l = px.^(-0.5);

deriv = -0.5*(px.^(-1)).*(c2*k2*(sech(k2*(x_bar - 1)).^2)/xd - (k1*c1*(1 + x_bar).^(-k1-1))/xd);

p = zeros(length(y),length(px));
bx = zeros(length(y),length(px));
by = zeros(length(y),length(px));
bz = zeros(length(y),length(px));
for i = 1:length(y)
    p(i,:) = px.*cosh(y(i)*l.^(-1)).^(-2) + 0.25;
    bx(i,:) = -l.^(-1).*tanh(y(i)*l.^(-1));
    by(i,:) = deriv.*(1 - y(i)*l.^(-1).*tanh(y(i)*l.^(-1)));
end

figure
starty = -10:10;
startx = 150*ones(length(starty),1);
[sx,sy] = meshgrid(x,y);
streamline(stream2(sx,sy,bx,by,startx,starty));
streamline(stream2(sx,sy,-bx,-by,startx,starty));
streamline(stream2(sx,sy,bx,by,[0 0 0],[-1.4 -1 -0.5]))


figure 
pcolor(p)