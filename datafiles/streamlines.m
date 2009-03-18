load solution.out
load nodes.out

fid=fopen('../domain');
a=textscan(fid,'%n%n%n%n','Delimiter',',');

up=a{1}(1);
down=a{2}(1);
left=a{3}(1);
right=a{4}(1);


fid=fopen('../holefile');
a=textscan(fid,'%2d');
nholes=a{1}(1);
frewind(fid)

for i=1:nholes
    a=textscan(fid,'%n%n%n%n%c','Headerlines',1,'Delimiter',',');
    h(i)=a{1}(1);
    k(i)=a{2}(1);
    r(i)=a{3}(1);
end
fclose(fid);

u=solution(:,1);
v=solution(:,2);
p=solution(:,3);

mag=sqrt((u.^2+v.^2));

x=nodes(:,1);
y=nodes(:,2);

[xi yi]=meshgrid(linspace(left,right),linspace(down,up));

magi = griddata(x,y,mag,xi,yi);
ui = griddata(x,y,u,xi,yi);
vi = griddata(x,y,v,xi,yi);
pint = griddata(x(p~=0),y(p~=0),p(p~=0),xi,yi);

%ui(sqrt((xi-.5).^2+(yi-.5).^2)<1)=0;
%vi(sqrt((xi-.5).^2+(yi-.5).^2)<1)=0;
%contour(xi(5:95,5:95),yi(5:95,5:95),pint(5:95,5:95),20)
%quiverc2wcmap(xi,yi,-ui,-vi,.3);
%hold on
%contour(xi,yi,magi,20)
contour(xi,yi,pint,20)
hold on
colorbar('NorthOutside')
%brighten(.25)
axis equal
axis([left right down up])
%colorbar('NorthOutside')

h1=streamslice(xi,yi,ui,vi,1.5);
set(h1,'color','k');


t=linspace(0,2*pi);

for i=1:nholes
    
    xf=r(i).*cos(t);
    yf=r(i).*sin(t);
 
    fill(xf+h(i),yf+k(i),'w')
end


%colormap([0 0 0])
%contourf(xi,yi,pint)

 
%  figure()
%  hold on 
%  t=linspace(0,2*pi);
%     
%  for i=1:nholes
%     
%     xf=r(i).*cos(t);
%     yf=r(i).*sin(t);
%  
%     fill(xf+h(i),yf+k(i),'w')
%  end
%  daspect([1 1 1]); view(2)
%  [verts averts] = streamslice(xi,yi,ui,vi,1.5); 
%  sl = streamline([verts averts]);
%  axis tight off;
%  set(sl,'Visible','off')
%  iverts = interpstreamspeed(xi,yi,ui,vi,verts,.05);
%  set(gca,'DrawMode','fast','Position',[0 0 1 1],'ZLim',[1 5.1])
%  set(gcf,'Color','black')
%  h=line;
%  streamparticles(h,iverts, 200, ...
%      'Animate',100,'FrameRate',40, ...
%      'MarkerSize',4,'MarkerFaceColor','yellow')