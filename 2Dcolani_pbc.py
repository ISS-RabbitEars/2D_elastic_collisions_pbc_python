import numpy as np
import math
import random as rnd
import matplotlib.pyplot as plt
from matplotlib import animation

class vec:
	def __init__(self,a,b):
		self.x=a
		self.y=b
	def mag(self):
		return math.sqrt(self.x**2+self.y**2)
	def hat(self):
		return vec(self.x/self.mag(), self.y/self.mag())
	def dot(self,r1):
		return self.x*r1.x+self.y*r1.y
	def add(self,r1):
		return vec(self.x+r1.x, self.y+r1.y)
	def __add__(self,r1):
		return self.add(r1)
	def sub(self,r1):
		return vec(self.x-r1.x, self.y-r1.y)
	def __sub__(self,r1):
		return self.sub(r1)
	def mul(self,a):
		return vec(a*self.x, a*self.y)
	def __mul__(self,a):
		return self.mul(a)
	def __rmul__(self,a):
		return self*a
	def div(self,a):
		return vec(self.x/a, self.y/a)
	def rot(self,a):
		return vec(self.x*math.cos(a)+self.y*math.sin(a), -self.x*math.sin(a)+self.y*math.cos(a))
	def __iadd__(self,r1):
		return self+r1
	def __isub__(self,r1):
		return self-r1

def cb(r1,rs1,s1,v1,t,n):
	for i in range(n):
		if r1[i].x<-s1:
			r1[i].x+=2*s1
		if r1[i].x>s1:
			r1[i].x-=2*s1
		if r1[i].y<-s1:
			r1[i].y+=2*s1
		if r1[i].y>s1:
			r1[i].y-=2*s1
	return r1,v1

def cc(r1,rs1,s1,v1,m1,t,n,col1):
	for i in range(n):
		drc=[]
		drj=[]
		drlmi=[]
		col=False
		for j in range(n):
			if j!=i:
				rt=r1[j]-r1[i]
				rtm=rt.mag()
				rss=rs1[i]+rs1[j]
				dx=r1[j].x-r1[i].x
				dy=r1[j].y-r1[i].y
				dxp1=r1[j].x-r1[i].x+2*s1
				dxp2=r1[j].x-r1[i].x-2*s1
				dyp1=r1[j].y-r1[i].y+2*s1
				dyp2=r1[j].y-r1[i].y-2*s1
				drxp1=math.sqrt(dxp1**2+dy**2)
				drxp2=math.sqrt(dxp2**2+dy**2)
				dryp1=math.sqrt(dx**2+dyp1**2)
				dryp2=math.sqrt(dx**2+dyp2**2)
				drxp1yp1=math.sqrt(dxp1**2+dyp1**2)
				drxp1yp2=math.sqrt(dxp1**2+dyp2**2)
				drxp2yp1=math.sqrt(dxp2**2+dyp1**2)
				drxp2yp2=math.sqrt(dxp2**2+dyp2**2)
				drl=[rtm,drxp1,drxp2,dryp1,dryp2,drxp1yp1,drxp1yp2,drxp2yp1,drxp2yp2]
				drlm=min(drl)
				if drlm<=rss:
					drc.append(drlm)
					drj.append(j)
					ii=drl.index(drlm)
					drlmi.append(ii)
					col=True
					col1=True
		if col:
			drmin=min(drc)
			ii=drc.index(drmin)
			j=drj[ii]
			id=drlmi[ii]
			r1[i]-=v1[i]*t
			r1[j]-=v1[j]*t
			if id==0:
				drx=r1[j].x-r1[i].x
				dry=r1[j].y-r1[i].y
			elif id==1:
				drx=r1[j].x-r1[i].x+2*s1
				dry=r1[j].y-r1[i].y
			elif id==2:
				drx=r1[j].x-r1[i].x-2*s1
				dry=r1[j].y-r1[i].y
			elif id==3:
				drx=r1[j].x-r1[i].x
				dry=r1[j].y-r1[i].y+2*s1
			elif id==4:
				drx=r1[j].x-r1[i].x
				dry=r1[j].y-r1[i].y-2*s1
			elif id==5:
				drx=r1[j].x-r1[i].x+2*s1
				dry=r1[j].y-r1[i].y+2*s1
			elif id==6:
				drx=r1[j].x-r1[i].x+2*s1
				dry=r1[j].y-r1[i].y-2*s1
			elif id==7:
				drx=r1[j].x-r1[i].x-2*s1
				dry=r1[j].y-r1[i].y+2*s1
			elif id==8:
				drx=r1[j].x-r1[i].x-2*s1
				dry=r1[j].y-r1[i].y-2*s1
			dvx=v1[j].x-v1[i].x
			dvy=v1[j].y-v1[i].y
			qa=dvx**2+dvy**2
			qb=2*(drx*dvx+dry*dvy)
			sr=rs1[i]+rs1[j]
			qc=drx**2+dry**2-sr**2
			dt1=(-qb-math.sqrt(qb**2-4*qa*qc))/(2*qa)
			dt2=t-dt1
			r1[i]+=v1[i]*dt1
			r1[j]+=v1[j]*dt1
			dr=vec(0,0)
			if id==0:
				dr.x=r1[j].x-r1[i].x
				dr.y=r1[j].y-r1[i].y
			elif id==1:
				dr.x=r1[j].x-r1[i].x+2*s1
				dr.y=r1[j].y-r1[i].y
			elif id==2:
				dr.x=r1[j].x-r1[i].x-2*s1
				dr.y=r1[j].y-r1[i].y
			elif id==3:
				dr.x=r1[j].x-r1[i].x
				dr.y=r1[j].y-r1[i].y+2*s1
			elif id==4:
				dr.x=r1[j].x-r1[i].x
				dr.y=r1[j].y-r1[i].y-2*s1
			elif id==5:
				dr.x=r1[j].x-r1[i].x+2*s1
				dr.y=r1[j].y-r1[i].y+2*s1
			elif id==6:
				dr.x=r1[j].x-r1[i].x+2*s1
				dr.y=r1[j].y-r1[i].y-2*s1
			elif id==7:
				dr.x=r1[j].x-r1[i].x-2*s1
				dr.y=r1[j].y-r1[i].y+2*s1
			elif id==8:
				dr.x=r1[j].x-r1[i].x-2*s1
				dr.y=r1[j].y-r1[i].y-2*s1
			drh=dr.hat()
			tr=math.acos(drh.x)
			if drh.y<0:
				tr=-tr
			vri=v1[i].rot(tr)
			vrj=v1[j].rot(tr)
			sm=m1[i]+m1[j]
			dm=m1[j]-m1[i]
			vifx=(-dm/sm)*vri.x+(2*m1[j]/sm)*vrj.x
			vjfx=(2*m1[i]/sm)*vri.x+(dm/sm)*vrj.x
			vri.x=vifx
			vrj.x=vjfx
			v1[i]=vri.rot(-tr)
			v1[j]=vrj.rot(-tr)
			r1[i]+=v1[i]*dt2
			r1[j]+=v1[j]*dt2
	if r1[i].x<-s1:
		r1[i].x+=2*s1
	if r1[i].y<-s1:
		r1[i].y+=2*s1
	if r1[i].x>s1:
		r1[i].x-=2*s1
	if r1[i].y>s1:
		r1[i].y-=2*s1
	if r1[j].x<-s1:
		r1[j].x+=2*s1
	if r1[j].y<-s1:
		r1[j].y+=2*s1
	if r1[j].x>s1:
		r1[j].x-=2*s1
	if r1[j].y>s1:
		r1[j].y-=2*s1
	return r1,v1

N=150
s=10
ma=1
mb=5
va=1
vb=5
dt=0.01
malpha=mb-ma
if malpha!=0:
	mbeta=ma/malpha
else:
	mbeta=0
valpha=vb-va
vbeta=va/valpha
r=[]
v=[]
m=[]
rs=[]

fig, ax = plt.subplots()

for i in range(N):
	if malpha!=0:
		m.append(malpha*(mbeta+rnd.random()))
	else:
		m.append(ma)
	t=2*np.pi*rnd.random()
	vm=valpha*(vbeta+rnd.random())
	v.append(vec(vm*math.cos(t),vm*math.sin(t)))
	rs.append(m[i]/(2*mb))
rx=(s-rs[0])*(1-2*rnd.random())
ry=(s-rs[0])*(1-2*rnd.random())
r.append(vec(rx,ry))
for i in range(N-1):
	col=True
	while col:
		col=False
		xt=(s-2*rs[i])*(1-2*rnd.random())
		yt=(s-2*rs[i])*(1-2*rnd.random())
		rt=vec(xt,yt)
		j=0
		while j<i+1 and not col:
			dr=rt-r[j]
			if dr.mag()<=(rs[i]+rs[j]):
				col=True
			j+=1
	r.append(rt)

def run(frame):
	global r,v
	for i in range(N):
		r[i]+=v[i]*dt
	pcol=True
	while pcol:
		pcol=False
		r,v=cb(r,rs,s,v,dt,N)
		r,v=cc(r,rs,s,v,m,dt,N,pcol)
	plt.clf()
	ax=plt.gca()
	ax.set_aspect(1)
	plt.xlim([-s,s])
	plt.ylim([-s,s])
	ax.set_facecolor('xkcd:black')
	for i in range(N):
		circle=plt.Circle((r[i].x,r[i].y),radius=rs[i],fc='r')
		plt.gca().add_patch(circle)
		if r[i].x-rs[i]<-s:
			circle = plt.Circle((r[i].x+2*s,r[i].y), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].x+rs[i]>s:
			circle = plt.Circle((r[i].x-2*s,r[i].y), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].y-rs[i]<-s:
			circle = plt.Circle((r[i].x,r[i].y+2*s), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].y+rs[i]>s:
			circle = plt.Circle((r[i].x,r[i].y-2*s), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].x-rs[i]<-s and r[i].y-rs[i]<-s:
			circle = plt.Circle((r[i].x+2*s,r[i].y+2*s), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].x-rs[i]<-s and r[s].y+rs[i]>s:
			circle = plt.Circle((r[i].x+2*s,r[i].y-2*s), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].x+rs[i]>s and r[i].y-rs[i]<-s:
			circle = plt.Circle((r[i].x-2*s,r[i].y+2*s), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)
		if r[i].x+rs[i]>s and r[i].y+rs[i]>s:
			circle = plt.Circle((r[i].x-2*s,r[i].y-2*s), radius=rs[i], fc='xkcd:red')
			plt.gca().add_patch(circle)


frps=60
sec=60
ani=animation.FuncAnimation(fig,run,frames=frps*sec)
#writervideo=animation.FFMpegWriter(fps=frps)
#ani.save('2dcol_pbc.mp4',writer=writervideo)
plt.show()


