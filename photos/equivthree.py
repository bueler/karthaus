#!/usr/bin/env python

#  make diagram which says "any two imply third" for geometry-change
#  equations

from pylab import *

def mycircle(R,x0,y0):
    theta = linspace(0.0,2.0*pi,400)
    x = x0 + R * cos(theta)
    y = y0 + R * sin(theta)
    return x,y

R = 0.7
FS = 24.0

x1,y1 = mycircle(R,-1.0,1.0)
plot(x1,y1,'k')
text(-1.0-0.5*R,1.0+0.2*R,'surface',fontsize=FS)
text(-1.0-0.8*R,1.0-0.2*R,'kinematical',fontsize=FS)

x2,y2 = mycircle(R,1.0,1.0)
plot(x2,y2,'k')
text(1.0-0.35*R,1.0+0.2*R,'base',fontsize=FS)
text(1.0-0.8*R,1.0-0.2*R,'kinematical',fontsize=FS)

x3,y3 = mycircle(R,0.0,-0.7)
plot(x3,y3,'k')
text(0.0-0.4*R,-0.7+0.2*R,'mass',fontsize=FS)
text(0.0-0.7*R,-0.7-0.2*R,'continuity',fontsize=FS)

axis('equal')
axis('off')

#show()
savefig('equivthree.pdf')

