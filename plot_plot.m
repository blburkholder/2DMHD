%run test first to get B

%filename = 'test.gif';
figure(1)
for where = 1:51

pcolor(x,y,reshape(p(where,:),[w,l])')
shading interp
colorbar

startx = 85:20:145;
starty = zeros(length(startx),1);
streamline(stream2(sx,sy,reshape(bx(where,:),[w,l])',reshape(by(where,:),[w,l])',startx,starty));
streamline(stream2(sx,sy,-reshape(bx(where,:),[w,l])',-reshape(by(where,:),[w,l])',startx,starty));

starty = 3:ddy/2;
startx = ddx*ones(length(starty),1);
[sx,sy] = meshgrid(x,y);
streamline(stream2(sx,sy,reshape(bx(where,:),[w,l])',reshape(by(where,:),[w,l])',startx,starty));
streamline(stream2(sx,sy,-reshape(bx(where,:),[w,l])',-reshape(by(where,:),[w,l])',startx,starty));

starty = 0:0.5:3;
startx = zeros(length(starty),1);
[sx,sy] = meshgrid(x,y);
streamline(stream2(sx,sy,reshape(bx(where,:),[w,l])',reshape(by(where,:),[w,l])',startx,starty));
streamline(stream2(sx,sy,-reshape(bx(where,:),[w,l])',-reshape(by(where,:),[w,l])',startx,starty));

%starty = 1:10;
%startx = zeros(length(starty),1);
%streamline(stream2(sx,sy,-reshape(bx(where,:),[l,w]),-reshape(by(where,:),[l,w]),startx,starty));
%starty = -10:-1;
%startx = zeros(length(starty),1);
%streamline(stream2(sx,sy,reshape(bx(where,:),[l,w]),reshape(by(where,:),[l,w]),startx,starty));

%starty = 2:5;
%startx = 150*ones(length(starty),1);
%streamline(stream2(sx,sy,reshape(bx(where,:),[l,w]),reshape(by(where,:),[l,w]),startx,starty));
%starty = -5:-2;
%startx = zeros(length(starty),1);
%streamline(stream2(sx,sy,-reshape(bx(where,:),[l,w]),-reshape(by(where,:),[l,w]),startx,starty));


startx = 5:10:75;
starty = zeros(length(startx),1);
for j = startx
    this = stream2(sx,sy,reshape(bx(where,:),[w,l])',reshape(by(where,:),[w,l])',j,0,[0.01,5000*(j-startx(1)+1)]);
    this_this = stream2(sx,sy,reshape(-bx(where,:),[w,l])',reshape(-by(where,:),[w,l])',j,0,[0.01,5000*(j-startx(1)+1)]);
    streamline(this);
    streamline(this_this);
end
pause(0.1)
%{
drawnow
frame = getframe(1)
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if where == 1
    imwrite(imind,cm,filename,'gif','Loopcount',inf);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append');
end
%}


end

hold on
res = 20;
uxx = reshape(ux(where,:),[l,w]);
uyy = reshape(uy(where,:),[l,w]); 
%quiver(x(1:res:end),y(1:res:end),uxx(1:res:end,1:res:end),uyy(1:res:end,1:res:end),'color','w')

title('t = ??')
xlabel('x')
ylabel('y')

