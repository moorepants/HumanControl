% example file take from http://bicycle.tudelft.nl/schwab/pshacker/ and
% modified slighty
global PSFID

PSFID=fopen(['plots' filesep 'benchmarkBicycle.ps'],'w')
header
boundingbox(152,169,405,357)
pstitle('Bicycle model degrees of freedom')
creator('pshacker')
creationdate(date)
prolog
dots=0.24;
translate(200,200)
selectfont('Helvetica-Oblique',10)

% scaling factor
scf=140;

% bicycle dimensions
l=scf*1.02;
Rr=scf*0.33;
Rf=scf*0.33;
lambda=pi/10;
trail=scf*0.08;
cmrfx=scf*0.3;
cmrfy=scf*0.9;
cmffx=scf*0.9;
cmffy=scf*0.7;

setlinewidth(8*dots)
setlinejoin(1)
setlinecap(1)

lwrf=12;
grrf=0.5;
lwff=12;
grff=0.75;

% newpath
% circle(0,Rr,Rr)
% circle(l,Rf,Rf)
% stroke

newpath
a0=39;
moveto(0+Rr*cosd(a0),Rr+Rr*sind(a0))
arc(0,Rr,Rr,a0,a0-12)
b0=119;
moveto(l+Rf*cosd(b0),Rf+Rf*sind(b0))
arc(l,Rf,Rf,b0,b0-12)
stroke

setlinewidth(lwrf*dots)
setgray(grrf)
newpath
moveto(0,Rr)
hh=scf*0.85;
py=hh;
px=l+trail-py*tan(lambda);
lineto(px,py)
moveto(cmrfx,cmrfx*(py-Rr)/px+Rr)
%moveto(cmrfx,39)
lineto(cmrfx,cmrfy)
stroke

setlinewidth(lwff*dots)
translate(px,py)
lambdad=lambda*180/pi;
rotate(lambdad)
%	bx=2;by=5;
%	bx=4;by=8;
	bx=4;by=7;
	setgray(1)
	rectfill(-bx,-by,2*bx,2*by)
	setgray(0)
    
    setlinewidth(lwrf*dots)
	setgray(grrf)
    newpath
	moveto(-bx,-by)
	lineto(-bx,by)
	moveto(bx,-by)
	lineto(bx,by)
    stroke

    uc = -trail*cos(lambda)+Rf*sin(lambda);
    vc = -py/cos(lambda)+trail*sin(lambda)+Rf*cos(lambda)-2;

    setlinewidth(lwff*dots)
    setgray(grff)
    newpath
	ls=20;hs=28;
	moveto(-ls,hs)
	lineto(0,hs)
    curveto(0,-60,-1,-50,uc,vc)
%	lineto(0,-(py-Rf)/cos(lambda))
% add stops on the front fork
    sx=bx-lwff*dots/2;
    sy=by+3.4;
    moveto(-sx,-sy)
    lineto( sx,-sy)
     moveto(-sx,sy)
     lineto(sx,sy)
	stroke


rotate(-lambdad)
translate(-px,-py)

setgray(0)

setlinewidth(2*dots)
dx=9;dy=3;lp=30;
newpath
pijl(0,0,l+lp+16,0,dx,dy)
pijl(0,0,0,-lp,dx,dy)
stroke

newpath
circle(0,0,2)
circle(l,0,2)
fil

setlinewidth(1*dots)
% newpath
% moveto(-Rr,0)
% lineto(0,0)
% moveto(lp+3+lpd+3+3,0)
% lineto(l+Rf,0)
% stroke

% Rpr=Rr-5;
% Rpf=Rf+5;
newpath
% moveto(0,0)
% lineto(0,2*Rr-5)
moveto(l,-23)
lineto(l,2*Rf-5)

% dimension arrows w and t
arrowfac=0.7;
zw=-lp+13;
pijl(0,zw,l,zw,arrowfac*dx,arrowfac*dy)
tx=px+hh*tan(lambda);
bt=-20;
moveto(tx,0)
lineto(tx,bt-3)
pijl(l,bt,tx,bt,arrowfac*dx,arrowfac*dy)
headangled=90-lambdad;
% rha=16;
% x1=tx-rha;y1=0;
% x2=tx-rha*cosd(headangled);y2=rha*sind(headangled);
rlam=Rf-5;
x1=l; y1=trail/tan(lambda);
x2=l-rlam*sin(lambda); y2=y1+rlam*cos(lambda);
moveto(x1,y1)
arc(x1,y1,rlam,90,90+lambdad)
pijlpunt(x2,y2,180+lambdad-2,arrowfac*dx,arrowfac*dy)
% moveto(x1,y1)
% arcn(tx,0,rha,180,180-headangled)
% pijlpunt(x1,y1,270-10,dx,dy)
% pijlpunt(x2,y2,90-headangled+10,dx,dy)
stroke

setgray(0)
translate(tx,0)
rotate(lambdad)
    setlinewidth(1*dots)
    newpath
    moveto(0,-18)
    lineto(0,3*Rf)
%     moveto(0,hs+14)
%     lineto(0,-hh/cos(lambda)-15)
    stroke
rotate(-lambdad)
translate(-tx,0)

rcm=6;
mcmmarker(0,Rr,rcm)
mcmmarker(l,Rf,rcm)
mcmmarker(cmrfx,cmrfy,rcm,grrf)
mcmmarker(cmffx,cmffy,rcm,grff)


% the symbols
translate(-200,-200)
selectfont('Helvetica-Oblique',10)
% setgray(1)
% rectfill(194,202,8,8)
setgray(0)
moveto(372,191)
show('x')
moveto(191,170)
show('z')
moveto(268,185)
show('w')
moveto(345,184)
show('c')
% setgray(1)
% rectfill(245,208,15,8)
% setgray(0)
selectfont('Symbol',10)
moveto(333,280)
show('l')
%selectfont('Helvetica-Oblique',8)
%rmoveto(0,-2)
%show('s')
selectfont('Helvetica',10)
moveto(192,190)
show('P')
moveto(333,190)
show('Q')
selectfont('Helvetica-Bold',10)
moveto(157,346)
show('A')

% body parts with names and symbols
selectfont('Helvetica',10)
moveto(158,297)
show('Rear wheel, R')
moveto(340,297)
show('Front wheel, F')
xB=176;
moveto(xB,335)
show('rear frame including')
cH=12;
moveto(xB,335-cH)
show('rider Body, B')
xH1=288;
moveto(xH1,349)
show('front frame \(fork and')
xH2=317;
moveto(xH2,349-cH)
show('Handlebar\), H')
moveto(361,178)
show('steer axis')



% selectfont('Symbol',8)
% moveto(175,150)
% show('abcdefghijklmnopqrstuvwxyz')
% selectfont('Helvetica-Oblique',8)
% moveto(175,140)
% show('abcdefghijklmnopqrstuvwxyz')

showpage
trailer
eof
fclose(PSFID)
