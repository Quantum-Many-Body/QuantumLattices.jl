import numpy as np
import matplotlib.pyplot as plt
import pdb

plt.ion()
plt.axis('off')
plt.axis('equal')

x1,x2,x3=0.0,1.0,0.5
y1,y2,y3=0.0,0.0,np.sqrt(3)/2
plt.plot([x1,x2],[y1,y2],lw=4,zorder=1)
plt.plot([x2,x3],[y2,y3],lw=4,zorder=1)
plt.plot([x3,x1],[y3,y1],lw=4,zorder=1)

c1=(89/255,90/255,15/255)
c2=(11/12,5/12,12/12)
c3=(12/(30+2/3),4/(30+2/3),27/(30+2/3))

plt.scatter(x1,y1,s=10000,color=c1,alpha=0.5)
plt.scatter(x2,y2,s=10000,color=c2,alpha=0.5)
plt.scatter(x3,y3,s=10000,color=c3,alpha=0.5)

r,theta=0.4,np.pi/9*3.5
dx=r*np.cos(theta)
dy=r*np.sin(theta)
plt.arrow(-dx,-dy,2*dx,2*dy,width=0.05,color='black',length_includes_head=True,zorder=3)
plt.arrow(1+dx,dy,-2*dx,-2*dy,width=0.05,color='black',length_includes_head=True,zorder=3)
plt.text(0.5,np.sqrt(3)/2-r*0.15,"?",fontsize=65,ha='center',va='center',zorder=3)
plt.xlim(-1.0,2.2)
plt.ylim(-1.5,2.2)

pdb.set_trace()
plt.savefig("log.pdf")
plt.close()
