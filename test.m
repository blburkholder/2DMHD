
bx = textread('bx.txt');
1
by = textread('by.txt');
2
%bz = textread('bz.txt');
p = textread('p.txt');
3
%p = textread('rho.txt');
ux = textread('ux.txt');
4
uy = textread('uy.txt');
5
%uz = textread('uz.txt');

bx = bx(:,1:end-1);
by = by(:,1:end-1);
%bz = bz(:,1:end-1);
p = p(:,1:end-1);
%rho = rho(:,1:end-1);
ux = ux(:,1:end-1);
uy = uy(:,1:end-1);
%uz = uz(:,1:end-1);
%}
%recon
dx = 0.5;
dy = 0.1;
ddx = 150;
ddy = 20;

%alfven
%dx = 0.5;
%dy = 0.5;
%ddx = 30;
%ddy = 30;

l = ddy/dy+1;
w = ddx/dx+1;

y = -ddy/2:dy:ddy/2;
x= 0:dx:ddx; 

%for recon
starty = -ddy/2:ddy/2;
%starty = -2:0.1:2;
startx = ddx*ones(length(starty),1);
[sx,sy] = meshgrid(x,y);

%for alven wave
%starty = -10:10;
%startx = zeros(length(starty),1);
%[sx,sy] = meshgrid(x,y);

[row,col] = size(p);

%movie

figure
for i = 1:row   
  %bxx = reshape(bx(i,:),[w,l])';
  %byy = reshape(by(i,:),[w,l])';  
  uxx = reshape(ux(i,:),[w,l])';
  uyy = reshape(uy(i,:),[w,l])';
  ppp = reshape(p(i,:),[w,l])';
 i
  pcolor(x,y,ppp)
  colorbar
  shading interp
  hold on
  quiver(x,y,uxx,uyy,'color','k')
  hold off
  
  streamline(stream2(sx,sy,reshape(bx(i,:),[w,l])',reshape(by(i,:),[w,l])',startx,starty));
  streamline(stream2(sx,sy,-reshape(bx(i,:),[w,l])',-reshape(by(i,:),[w,l])',startx,starty));
  
  pause
  %clf
end
%}

%{
figure
pcolor(x,y,reshape(p(1,:),[l,w]))
shading interp
colorbar
starty = -ddy/2:ddy/2;
startx = ddx*ones(length(starty),1);
streamline(stream2(sx,sy,reshape(bx(1,:),[l,w]),reshape(by(1,:),[l,w]),startx,starty));
streamline(stream2(sx,sy,-reshape(bx(1,:),[l,w]),-reshape(by(1,:),[l,w]),startx,starty));
streamline(stream2(sx,sy,reshape(bx(1,:),[l,w]),reshape(by(1,:),[l,w]),[0 0],[-1.4 -1]));
streamline(stream2(sx,sy,reshape(bx(1,:),[l,w]),reshape(by(1,:),[l,w]),[ddx],[0.6]));

starty = 6:10;
startx = zeros(length(starty),1);
streamline(stream2(sx,sy,-reshape(bx(1,:),[l,w]),-reshape(by(1,:),[l,w]),startx,starty));

starty = -10:-6;
streamline(stream2(sx,sy, reshape(bx(1,:),[l,w]),reshape(by(1,:),[l,w]),startx,starty));

title('initial condition')
xlabel('x')
ylabel('y')

plot_plot();
%}

