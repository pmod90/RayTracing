
#include <stdlib.h>

#include <windows.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <glut.h>
#endif
#include <conio.h>
#include<iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
using namespace std;

struct Point
{
   float x, y, z;
   float Kar,Kag,Kab;
   float Kdr,Kdg,Kdb;
   float Ksr,Ksg,Ksb;
   float Iar,Iag,Iab;
   float Isr,Isg,Isb;
   float n,flag;
   float spec,ldotn,tval,distance;
};


Point p,n,ls,l,temp1,R,V,sphere1,sphere2,cylinder,sphere11,sphere22,final1,plane,cube,cone,cyltop,cubeedge1,cubeedge2,cubeedge3,cubeedge4,cubeedge5,cubeedge6,p0,p1,p2,p3,p4,p10,p20,p30,p40,p12,p23,p34,p41;



float a,b,c,t0,t1,k,t,Ir,Ig,Ib;
float X0 = 0;
float Y0 = 0;
float Z0 = 0;
float xp;
float yp;
float RV;
float zp = -40;
float xleft = -40;
float xright = 40;
float ybottom = -40;
float ytop = 40;
float len,i,j;
float dx = 0.1;
float dy = 0.1;
float Xd,Yd,Zd;
std::vector<float> ut;
Point cross(Point a,Point b)
{
	//(a2b3−a3b2)i−(a1b3−a3b1)j+(a1b2−a2b1)k.
	Point c;
	c.x = (a.y*b.z) - (a.z*b.y);
	c.y = (-a.x*b.z) + (a.z*b.x);
	c.z = (a.x*b.y) - (a.y*b.x);
	return(c);
}
Point giveVector(Point a,Point b)
{
	Point c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return(c);	
}
float mag(float x1,float y1,float z1)
{float k1 = sqrt(x1*x1 + y1*y1 + z1*z1); 
	return(k1);
}
Point normalize(Point a)
{
	float mg = mag(a.x,a.y,a.z);
	a.x = a.x/mg;
	a.y = a.y/mg;
	a.z = a.z/mg;	
	return(a);	
}
float dot(float x11,float y11,float z11,float x22,float y22,float z22)
{float d1 = (x11*x22 + y11*y22 + z11*z22); 
	return(d1);
}
Point planeintersect(float Xc,float Yc,float Zc,int i,int j,Point plane,int orientation,float h)
{	
	float mg;
	/*Orientations
	1.XZ plane
	2.XY plane
	3.YZ plane
	*/
	if(orientation == 1)//xz
		{
			
	 p0.x = X0;
	 p0.y = Y0;
	 p0.z = Z0;
     		
	 p1.x = Xc + h/2.0;
	 p1.y = Yc;
	 p1.z = Zc - h/2.0;

	 p2.x = Xc + h/2.0;
	 p2.y = Yc;
	 p2.z = Zc + h/2.0;

	 p3.x = Xc - h/2.0;
	 p3.y = Yc;
	 p3.z = Zc + h/2.0;

	 p4.x = Xc - h/2.0;
	 p4.y = Yc;
	 p4.z = Zc - h/2.0;

	 p10 = giveVector(p1,p0);
	 p20 = giveVector(p2,p0);
	 p30 = giveVector(p3,p0);
	 p40 = giveVector(p4,p0);

	 n = cross(giveVector(p2,p1),giveVector(p2,p3));
	 n = normalize(n);

	 p12 = cross(p10,p20);
	 p23 = cross(p20,p30);
	 p34 = cross(p30,p40);
	 p41 = cross(p40,p10);


	}
	else if(orientation == 2)
	{
	/* n.x = 0;
	 n.y = 0;
	 n.z = 1;*/

	 p0.x = X0;
	 p0.y = Y0;
	 p0.z = Z0;
     		
	 p1.x = Xc + h/2.0;
	 p1.y = Yc + h/2.0;
	 p1.z = Zc;

	 p2.x = Xc + h/2.0;
	 p2.y = Yc - h/2.0;
	 p2.z = Zc;

	 p3.x = Xc - h/2.0;
	 p3.y = Yc - h/2.0;
	 p3.z = Zc;

	 p4.x = Xc - h/2.0;
	 p4.y = Yc + h/2.0;
	 p4.z = Zc;

	 p10 = giveVector(p1,p0);
	 p20 = giveVector(p2,p0);
	 p30 = giveVector(p3,p0);
	 p40 = giveVector(p4,p0);

	 p12 = cross(p10,p20);
	 p23 = cross(p20,p30);
	 p34 = cross(p30,p40);
	 p41 = cross(p40,p10);

	 	 n = cross(giveVector(p2,p1),giveVector(p2,p3));
	 n = normalize(n);


	}
	else if(orientation == 3)
	{
	/*n.x = 1;
	 n.y = 0;
	 n.z = 0;*/

	 p0.x = X0;
	 p0.y = Y0;
	 p0.z = Z0;
     		
	 p1.x = Xc;
	 p1.y = Yc + h/2.0;
	 p1.z = Zc - h/2.0;

	 p2.x = Xc;
	 p2.y = Yc - h/2.0;
	 p2.z = Zc - h/2.0;

	 p3.x = Xc;
	 p3.y = Yc - h/2.0;
	 p3.z = Zc + h/2.0;;

	 p4.x = Xc;
	 p4.y = Yc + h/2.0;
	 p4.z = Zc + h/2.0;

	 p10 = giveVector(p1,p0);
	 p20 = giveVector(p2,p0);
	 p30 = giveVector(p3,p0);
	 p40 = giveVector(p4,p0);


	 p12 = cross(p10,p20);
	 p23 = cross(p20,p30);
	 p34 = cross(p30,p40);
	 p41 = cross(p40,p10);
	
	 	 n = cross(giveVector(p2,p1),giveVector(p2,p3));
	 n = normalize(n);
	}
	xp = xleft + 0.5*dx + i*dx;
		yp = ybottom + 0.5*dy + j*dy - 20;
	
	Xd = xp - X0;
	Yd = yp - Y0;
	Zd = zp - Z0;
	len = sqrt(Xd*Xd + Yd*Yd + Zd*Zd);
	Xd = Xd/len;
	Yd = Yd/len;
	Zd = Zd/len;
	float d = -n.x*Xc -n.y*Yc - n.z*Zc;
	if(d>0.0)
	{t = -(n.x*X0 + n.y*Y0 + n.z*Z0 + d)/(n.x*Xd + n.y*Yd + n.z*Zd);
	if(t>0.0)
	{
	

	 float d1 = dot(p12.x,p12.y,p12.z,Xd,Yd,Zd);
	 float d2 = dot(p23.x,p23.y,p23.z,Xd,Yd,Zd);
	 float d3 = dot(p34.x,p34.y,p34.z,Xd,Yd,Zd);
	 float d4 = dot(p41.x,p41.y,p41.z,Xd,Yd,Zd);
	
	 if((d1>0.0)&&(d2>0.0)&&(d3>0.0)&&(d4>0.0))
	 {
		plane.x = X0 + Xd*t;
		plane.y = Y0 + Yd*t;
		plane.z = Z0 + Zd*t;
		ls.x = 20;
		ls.y = 20;
		ls.z = -5;
		l.x = ls.x - plane.x;
		l.y = ls.y - plane.y;
		l.z = ls.z - plane.z;
		 mg = mag(l.x,l.y,l.z);
		//            N = ((x - cx)/R, (y - cy)/R, (z - cz)/R)
		 l.x = l.x/mg;
		 l.y = l.y/mg;
		 l.z = l.z/mg;
		 plane.ldotn = dot(n.x,n.y,n.z,l.x,l.y,l.z);
		 temp1.x = 2*(n.x)*(plane.ldotn);
		 temp1.y = 2*(n.y)*(plane.ldotn);
		 temp1.z = 2*(n.z)*(plane.ldotn);
		 R.x = temp1.x - l.x;
		 R.y = temp1.y - l.y; 
		 R.z = temp1.z - l.z;
		  mg = mag(R.x,R.y,R.z);
		  R.x = R.x/mg;
		 R.y = R.y/mg;
		 R.z = R.z/mg;
		 V.x = X0 - plane.x;
		 V.y = Y0 - plane.y;
	     V.z = Z0 - plane.z;
		 mg = mag(V.x,V.y,V.z);
		 
		 V.x = V.x/mg;
		 V.y = V.y/mg;
		 V.z = V.z/mg;
		 RV = dot(R.x,R.y,R.z,V.x,V.y,V.z);
		  plane.spec = powf(RV,plane.n);
		  plane.flag = 0.0;
		  plane.tval = 1.0;
		  plane.distance = mag(plane.x,plane.y,plane.z);
		  if(h!=12)
		  {return plane;}
		  else if(((plane.x - Xc)*(plane.x - Xc) + (plane.y - Yc)*(plane.y - Yc) + (plane.z - Zc)*(plane.z - Zc) - (h/2.0)*(h/2.0)) <= 0.0)
		  {return plane;}
		  else
		  {
		  plane.x = 0.0;
	 plane.y = 0.0;
	 plane.z = 0.0;
	 plane.flag = 1.0;
	 plane.tval = 0.0;
	 plane.distance = 1000.0;
	 return plane;}
		  

	 }
	 else
	 {
		plane.x = 0.0;
	 plane.y = 0.0;
	 plane.z = 0.0;
	 plane.flag = 1.0;
	 plane.tval = 0.0;
	 plane.distance = 1000.0;
	 return plane;


	 }
	}else{
	
		plane.x = 0.0;
	 plane.y = 0.0;
	 plane.z = 0.0;
	 plane.flag = 1.0;
	 plane.tval = 0.0;
	 plane.distance = 1000.0;
	 return plane;
}}
	else
	{
		
		plane.x = 0.0;
	 plane.y = 0.0;
	 plane.z = 0.0;
	 plane.flag = 1.0;
	 plane.tval = 0.0;
	 plane.distance = 1000.0;
	 return plane;

	}}
Point coneintersect(float Xc,float Yc,float Zc,int i,int j,Point cone,float r,float h)
{
//Actual s;
//s.intersection = setpoint(0.0, 0.0, 0.0);
	float mg;
xp = xleft + 0.5*dx + i*dx;
		yp = ybottom + 0.5*dy + j*dy - 20;
		Xd = xp - X0;
		Yd = yp - Y0;
		Zd = zp - Z0;
		len = sqrt(Xd*Xd + Yd*Yd + Zd*Zd);
		Xd = Xd/len;
		Yd = Yd/len;
		Zd = Zd/len;
		if((j == 201)&&(i == 799))
		{int d = 2;}
	float r1 = r*r;
float h1 = h*h;
float rh = r1/h1;
//float t11 = Xd*(X0-Xc);
//float t22 = Zd*(Z0-Zc);
//float t32 = Yd*(Y0-Yc-h)*rh;
float a = Xd*Xd + Zd*Zd - Yd*Yd*rh;
float b = 2*(Xd*(X0-Xc) + Zd*(Z0-Zc) - Yd*(Y0-Yc-h)*rh)/a;
//float b = 2*(p.x*d.x-cr.x*d.x + p.z*d.z-cr.z*d.z-(p.y-cr.y-h)*d.y*rh);
float c = ((X0-Xc)*(X0-Xc) + (Z0-Zc)*(Z0-Zc) - ((Y0-Yc)*(Y0-Yc) - 2*h*(Y0-Yc))*rh - r*r)/a;
//float c = (p.x-cr.x)*(p.x-cr.x) + (p.z-cr.z)*(p.z-cr.z) - ((p.y-cr.y)*(p.y-cr.y)-2*h*p.y+ 2*cr.y*h)*rh-r*r;
float t = 0;
float t0 = 0;
float t1 = 0;
float k = b*b;
float dis = k - 4*c;

if(dis==0)
{

t=-b/2;
}
else if(dis>0)
{
	
t0 = (-b-sqrt(dis))/2;
t1 = (-b+sqrt(dis))/2;
if (t0>0 && t1==0)
{
t=t0;
}
else if(t1>0 && t0==0)
{
t=t1;
}
else if (t0>0 && t1>0)
{
if (t0>t1)
{
t=t1;
}
else
{
t=t0;
}
}
else if (t0>0 && t1<0)
{
t = t0;
}
else if (t0<0 && t1>0)
{
t = t1;
}
cone.x = (X0 + Xd*t);
		cone.y = (Y0 + Yd*t);
		cone.z = (Z0 + Zd*t);
		/*

		float rl = r/h;
float dx = i.x - c.x;
float dz = i.z - c.z;
float d = sqrt(dx*dx +dz*dz);
float t = i.y - d*tan(rl);
Point sn = setpoint((i.x - c.x), (i.y - t), (i.z -c.z));
return sn;
*/
float rl = r/h;

float dX = cone.x - Xc;
float dZ = cone.z - Zc;
float a = sqrt(dX*dX +dZ*dZ);
float b = cone.y - a*tan(r/h);
n.x = cone.x - Xc;
n.y = cone.y - b;
n.z = cone.z - Zc;
n = normalize(n);

ls.x = 20;
		ls.y = 20;
		ls.z = -5;
		l.x = ls.x - cone.x;
		l.y = ls.y - cone.y;
		l.z = ls.z - cone.z;
		 mg = mag(l.x,l.y,l.z);
		//            N = ((x - cx)/R, (y - cy)/R, (z - cz)/R)
		 l.x = l.x/mg;
		 l.y = l.y/mg;
		 l.z = l.z/mg;
		 cone.ldotn = dot(n.x,n.y,n.z,l.x,l.y,l.z);
		  temp1.x = 2*(n.x)*(cone.ldotn);
		 temp1.y = 2*(n.y)*(cone.ldotn);
		 temp1.z = 2*(n.z)*(cone.ldotn);
		 R.x = temp1.x - l.x;
		 R.y = temp1.y - l.y; 
		 R.z = temp1.z - l.z;
		  mg = mag(R.x,R.y,R.z);
		  R.x = R.x/mg;
		 R.y = R.y/mg;
		 R.z = R.z/mg;
		 V.x = X0 - cone.x;
		 V.y = Y0 - cone.y;
	     V.z = Z0 - cone.z;
		 mg = mag(V.x,V.y,V.z);
		 V.x = V.x/mg;
		 V.y = V.y/mg;
		 V.z = V.z/mg;
		 RV = dot(R.x,R.y,R.z,V.x,V.y,V.z);
		  cone.spec = powf(RV,cone.n);
		  cone.flag = 0.0;
		  cone.tval = 1.0;
		  cone.distance = mag(cone.x,cone.y,cone.z);
 if((cone.y <= Yc)||(cone.y > (Yc + 8.0)))
		  { cone.x = 0.0;
	 cone.y = 0.0;
	 cone.z = 0.0;
	 cone.flag = 1.0;
	 cone.tval = 0.0;
	 cone.distance = 1000;
	 return cone;}
		  else
		  {
			  
	return cone;
		 }
}
else if(dis<0)
{cone.x = 0.0;
	 cone.y = 0.0;
	 cone.z = 0.0;
	 cone.flag = 1.0;
	 cone.tval = 0.0;
	 cone.distance = 1000;
	 return cone;}

}
Point cylinderintersect(float r,float h,float Xc1,float Yc1,float Zc1,int i,int j,Point cylinder)
{float mg;
xp = xleft + 0.5*dx + i*dx;
		yp = ybottom + 0.5*dy + j*dy - 20;
		Xd = xp - X0;
		Yd = yp - Y0;
		Zd = zp - Z0;
		len = sqrt(Xd*Xd + Yd*Yd + Zd*Zd);
		Xd = Xd/len;
		Yd = Yd/len;
		Zd = Zd/len;
		a = Xd*Xd + Zd*Zd;

b = 2*(Xd*(X0 - Xc1) + Zd*(Z0 - Zc1));
c = (X0 - Xc1)*(X0 - Xc1) + (Z0 - Zc1)*(Z0 - Zc1) - r*r;
float k = b*b;
float dis = k - 4*a*c;
if(dis==0)
{

t=-b/2;
}
else if(dis>0)
{

float t0 = (-b-sqrtf(dis))/2;
float t1 = (-b+sqrtf(dis))/2;
if (t0>0 && t1>0)
{
if (t0>t1)
{
t=t1;
}
else
{
t=t0;
}
}
else if (t0>0 && t1<0)
{
t = t0;
}
else if (t0<0 && t1>0)
{
t = t1;
}


cylinder.x = (X0 + Xd*t);
		cylinder.y = (Y0 + Yd*t);
		cylinder.z = (Z0 + Zd*t);
		n.x = (cylinder.x - Xc1);
		n.y = (cylinder.y - Yc1);
		n.z = (cylinder.z - Zc1);
		 mg = mag(n.x,n.y,n.z);
		 n.x = (cylinder.x - Xc1)/r;
		n.y = (cylinder.y - Yc1)/r;
		n.z = (cylinder.z - Zc1)/r;
		ls.x = 20;
		ls.y = 20;
		ls.z = -5;
		l.x = ls.x - cylinder.x;
		l.y = ls.y - cylinder.y;
		l.z = ls.z - cylinder.z;
		 mg = mag(l.x,l.y,l.z);
		//            N = ((x - cx)/R, (y - cy)/R, (z - cz)/R)
		 l.x = l.x/mg;
		 l.y = l.y/mg;
		 l.z = l.z/mg;
		 cylinder.ldotn = dot(n.x,n.y,n.z,l.x,l.y,l.z);
		  temp1.x = 2*(n.x)*(cylinder.ldotn);
		 temp1.y = 2*(n.y)*(cylinder.ldotn);
		 temp1.z = 2*(n.z)*(cylinder.ldotn);
		 R.x = temp1.x - l.x;
		 R.y = temp1.y - l.y; 
		 R.z = temp1.z - l.z;
		  mg = mag(R.x,R.y,R.z);
		  R.x = R.x/mg;
		 R.y = R.y/mg;
		 R.z = R.z/mg;
		 V.x = X0 - cylinder.x;
		 V.y = Y0 - cylinder.y;
	     V.z = Z0 - cylinder.z;
		 mg = mag(V.x,V.y,V.z);
		 V.x = V.x/mg;
		 V.y = V.y/mg;
		 V.z = V.z/mg;
		 RV = dot(R.x,R.y,R.z,V.x,V.y,V.z);
		  cylinder.spec = powf(RV,cylinder.n);
		  cylinder.flag = 0.0;
		  cylinder.tval = 1.0;
		  cylinder.distance = mag(cylinder.x,cylinder.y,cylinder.z);
		 if(cylinder.y <= Yc1 +2.5 || cylinder.y >= Yc1 + 6)
		  {cylinder.x = 0.0;
	 cylinder.y = 0.0;
	 cylinder.z = 0.0;
	 cylinder.flag = 1.0;
	 cylinder.tval = 0.0;
	 cylinder.distance = 1000;
	 return cylinder;}
		  else
		  {
	 return cylinder;
		 }}
else if(dis<0)
{cylinder.x = 0.0;
	 cylinder.y = 0.0;
	 cylinder.z = 0.0;
	 cylinder.flag = 1.0;
	 cylinder.tval = 0.0;
	 cylinder.distance = 1000;
	 return cylinder;}

}
Point sphereintersectdraw(float r,float Xc1,float Yc1,float Zc1,Point sphere,int i,int j)
{float rd = r;

	
	{	
		{xp = xleft + 0.5*dx + i*dx;
		yp = ybottom + 0.5*dy + j*dy - 20;
	
		Xd = xp - X0;
		Yd = yp - Y0;
		Zd = zp - Z0;
		len = sqrt(Xd*Xd + Yd*Yd + Zd*Zd);
		Xd = Xd/len;
		Yd = Yd/len;
		Zd = Zd/len;
			a = Xd*Xd + Yd*Yd + Zd*Zd;
			
b = 2*(Xd*(X0 - Xc1) + Yd*(Y0 - Yc1) + Zd*(Z0 - Zc1));

c = (X0 - Xc1)*(X0 - Xc1) + (Y0 - Yc1)*(Y0 - Yc1) + (Z0 - Zc1)*(Z0 - Zc1) - rd*rd;
k = b*b;
float disc = k - 4*a*c;
if (disc > 0)
       {   float distSqrt = sqrtf(disc);
t0 = (-b - distSqrt)/2.0*a;
t1 = (-b + distSqrt)/2.0*a;
 if (t0 > t1)
    {
        // if t0 is bigger than t1 swap them around
        float temp = t0;
        t0 = t1;
        t1 = temp;
    }
	 if (t1 < 0)
	 {sphere.x = 0.0;
	 sphere.y = 0.0;
	 sphere.z = 0.0;
	 sphere.flag = 1.0;
	 sphere.tval = 0.0;
	 sphere.distance = 1000.0;
	 
	 return sphere;}
	  if (t0 < 0)
    {float mg;
        t = t1;
		sphere.x = (X0 + Xd*t);
		sphere.y = (Y0 + Yd*t);
		sphere.z = (Z0 + Zd*t);
		n.x = (sphere.x - Xc1);
		n.y = (sphere.y - Yc1);
		n.z = (sphere.z - Zc1);
		 mg = mag(n.x,n.y,n.z);
		 n.x = (sphere.x - Xc1)/rd;
		n.y = (sphere.y - Yc1)/rd;
		n.z = (sphere.z - Zc1)/rd;
		ls.x = 20;
		ls.y = 20;
		ls.z = -5;
		l.x = ls.x - sphere.x;
		l.y = ls.y - sphere.y;
		l.z = ls.z - sphere.z;
		 mg = mag(l.x,l.y,l.z);
		//            N = ((x - cx)/R, (y - cy)/R, (z - cz)/R)
		 l.x = l.x/mg;
		 l.y = l.y/mg;
		 l.z = l.z/mg;
		 sphere.ldotn = dot(n.x,n.y,n.z,l.x,l.y,l.z);
		 temp1.x = 2*(n.x)*(sphere.ldotn);
		 temp1.y = 2*(n.y)*(sphere.ldotn);
		 temp1.z = 2*(n.z)*(sphere.ldotn);
		 R.x = temp1.x - l.x;
		 R.y = temp1.y - l.y; 
		 R.z = temp1.z - l.z;
		  mg = mag(R.x,R.y,R.z);
		  R.x = R.x/mg;
		 R.y = R.y/mg;
		 R.z = R.z/mg;
		 V.x = X0 - sphere.x;
		 V.y = Y0 - sphere.y;
	     V.z = Z0 - sphere.z;
		 mg = mag(V.x,V.y,V.z);
		 V.x = V.x/mg;
		 V.y = V.y/mg;
		 V.z = V.z/mg;
		 RV = dot(R.x,R.y,R.z,V.x,V.y,V.z);
		  sphere.spec = powf(RV,sphere.n);
		  sphere.flag = 0.0;
		  sphere.tval = 1.0;
		  sphere.distance = mag(sphere.x,sphere.y,sphere.z);
		 return sphere;

       
    }
    // else the intersection point is at t0
    else
    {float mg;
        t = t0;
	    sphere.x = (X0 + Xd*t);
		sphere.y = (Y0 + Yd*t);
		sphere.z = (Z0 + Zd*t);
		n.x = (sphere.x - Xc1);
		n.y = (sphere.y - Yc1);
		n.z = (sphere.z - Zc1);
		 mg = mag(n.x,n.y,n.z);
		 n.x = (sphere.x - Xc1)/rd;
		n.y = (sphere.y - Yc1)/rd;
		n.z = (sphere.z - Zc1)/rd;
		ls.x = 20;
		ls.y = 20;
		ls.z = -5;
		l.x = ls.x - sphere.x;
		l.y = ls.y - sphere.y;
		l.z = ls.z - sphere.z;
		 mg = mag(l.x,l.y,l.z);
		//            N = ((x - cx)/R, (y - cy)/R, (z - cz)/R)
		 l.x = l.x/mg;
		 l.y = l.y/mg;
		 l.z = l.z/mg;
		 sphere.ldotn = dot(n.x,n.y,n.z,l.x,l.y,l.z);
		  temp1.x = 2*(n.x)*(sphere.ldotn);
		 temp1.y = 2*(n.y)*(sphere.ldotn);
		 temp1.z = 2*(n.z)*(sphere.ldotn);
		 R.x = temp1.x - l.x;
		 R.y = temp1.y - l.y; 
		 R.z = temp1.z - l.z;
		  mg = mag(R.x,R.y,R.z);
		  R.x = R.x/mg;
		 R.y = R.y/mg;
		 R.z = R.z/mg;
		 V.x = X0 - sphere.x;
		 V.y = Y0 - sphere.y;
	     V.z = Z0 - sphere.z;
		 mg = mag(V.x,V.y,V.z);
		 V.x = V.x/mg;
		 V.y = V.y/mg;
		 V.z = V.z/mg;
		 RV = dot(R.x,R.y,R.z,V.x,V.y,V.z);
		  sphere.spec = powf(RV,sphere.n);
		  sphere.flag = 0.0;
		  sphere.tval = 1.0;
		  sphere.distance = mag(sphere.x,sphere.y,sphere.z);

		  return sphere;		 

    }

}	
else
{
	 sphere.x = 0.0;
	 sphere.y = 0.0;
	 sphere.z = 0.0;
	 sphere.flag = 1.0;
	 sphere.tval = 0.0;
	 sphere.distance = 1000;
	 return sphere;
}
	}}
}
void trace()
{for(i = 0;i<800;i++)
{for(j = 0;j<800;j++)
{
	plane = planeintersect(0.0,-20.0,-60.0,i,j,plane,1,80);
	sphere22 = sphereintersectdraw(8.0,-8.0,-12.0,-70.0,sphere2,i,j);
	sphere11 = sphereintersectdraw(12.0,15.0,-5.0,-80.0,sphere1,i,j);
	cylinder = cylinderintersect(6.0,6.0,-20.0,-20.0,-50.0,i,j,cylinder);
	cubeedge1 = planeintersect(0.0,-20.0,-50.0,i,j,cube,1,10);
	cubeedge2 = planeintersect(0.0,-10.0,-50.0,i,j,cube,1,10);
	cubeedge3 = planeintersect(-5.0,-15.0,-50.0,i,j,cube,3,10);
	cubeedge4 = planeintersect(5.0,-15.0,-50.0,i,j,cube,3,10);
	cubeedge5 = planeintersect(0.0,-15.0,-45.0,i,j,cube,2,10);
	cyltop = planeintersect(-20.0,-15.5,-50.0,i,j,cylinder,1,12);
	cone = coneintersect(20.0,-20.0,-50.0,i,j,cone,8,8);
	/*cubeedge6 = planeintersect(0.0,-15.0,-35.0,i,j,cube,2,10);*/

	if((cone.tval == 1.0)||(sphere22.tval==1.0)||(sphere11.tval==1.0)||(cylinder.tval==1.0)||(plane.tval==1.0)||(cubeedge1.tval==1.0)||(cubeedge2.tval==1.0)||(cubeedge3.tval==1.0)||(cubeedge4.tval==1.0)||(cubeedge5.tval==1.0)||(cyltop.tval==1.0))
	{
		
		float finaldistance = min(min(min(min(min(min(min(min(min(min(sphere11.distance,sphere22.distance),cylinder.distance),plane.distance),cubeedge1.distance),cubeedge2.distance),cubeedge3.distance),cubeedge4.distance),cubeedge5.distance),cyltop.distance),cone.distance);
		if((finaldistance)==sphere11.distance)
		{final1 = sphere11;}
		else if((finaldistance)==sphere22.distance)
		{final1 = sphere22;}
		else if((finaldistance)==cylinder.distance)
		{final1 = cylinder;}
		else if((finaldistance)==plane.distance)
		{final1 = plane;}
		else if((finaldistance)==cubeedge1.distance)
		{final1 = cubeedge1;}
		else if((finaldistance)==cubeedge2.distance)
		{final1 = cubeedge2;}
		else if((finaldistance)==cubeedge3.distance)
		{final1 = cubeedge3;}
		else if((finaldistance)==cubeedge4.distance)
		{final1 = cubeedge4;}
		else if((finaldistance)==cubeedge5.distance)
		{final1 = cubeedge5;}
	else if((finaldistance)==cyltop.distance)
		{final1 = cyltop;}
	else if((finaldistance)==cone.distance)
	{final1 = cone;
	}
	
	

		 Ir = ((final1.Iar)*(final1.Kar)) + ((final1.Isr)*(final1.Kdr)*(final1.ldotn)) + ((final1.Ksr)*(final1.spec));
		 Ig = ((final1.Iag)*(final1.Kag)) + ((final1.Isg)*(final1.Kdg)*(final1.ldotn)) + ((final1.Ksg)*(final1.spec));
		 Ib = ((final1.Iab)*(final1.Kab)) + ((final1.Isb)*(final1.Kdb)*(final1.ldotn)) + ((final1.Ksb)*(final1.spec));
	glColor3f(Ir,Ig,Ib);
	//if(final1.flag!=1.0)
 { glBegin(GL_POINTS);
		  glVertex2f(i,j);
		  glEnd();}}


}}}
void display(void) {
	/* clear the screen to the clear colour */
	glClear(GL_COLOR_BUFFER_BIT);
	sphere1.Kar = 0.231;
	sphere1.Kag = 0.231;
	sphere1.Kab = 0.231;
	sphere1.Kdr = 0.278;
	sphere1.Kdg = 0.278;
	sphere1.Kdb = 0.278;
	sphere1.Ksr = 0.774;
	sphere1.Ksg = 0.774;
	sphere1.Ksb = 0.774;
	sphere1.n = 89.0;
	sphere1.Iar = 1.0;
	sphere1.Iag = 1.0;
	sphere1.Iab = 1.0;
	sphere1.Isr = 1.0;
	sphere1.Isg = 1.0;
	sphere1.Isb = 1.0;
	sphere1.flag = 1.0;

	sphere2.Kar = 0.3;
	sphere2.Kag = 0.16;
	sphere2.Kab = 0.12;
	sphere2.Kdr = 0.8;
	sphere2.Kdg = 0.4;
	sphere2.Kdb = 0.35;
	sphere2.Ksr = 0.1;
	sphere2.Ksg = 0.1;
	sphere2.Ksb = 0.1;
	sphere2.n = 5.0;
	sphere2.Iar = 1.0;
	sphere2.Iag = 1.0;
	sphere2.Iab = 1.0;
	sphere2.Isr = 1.0;
	sphere2.Isg = 1.0;
	sphere2.Isb = 1.0;
	sphere2.flag = 2.0;

	cylinder.Kar = 0.35;
	cylinder.Kag = 0.06;
	cylinder.Kab = 0.02;
	cylinder.Kdr = 0.82;
	cylinder.Kdg = 0.2;
	cylinder.Kdb = 0.1;
	cylinder.Ksr = 0.1;
	cylinder.Ksg = 0.1;
	cylinder.Ksb = 0.1;
	cylinder.n = 5.0;
	cylinder.Iar = 1.0;
	cylinder.Iag = 1.0;
	cylinder.Iab = 1.0;
	cylinder.Isr = 1.0;
	cylinder.Isg = 1.0;
	cylinder.Isb = 1.0;

	plane.Kar = 0.3;
	plane.Kag = 0.3;
	plane.Kab = 0.3;
	plane.Kdr = 0.6;
	plane.Kdg = 0.6;
	plane.Kdb = 0.6;
	plane.Ksr = 0.1;
	plane.Ksg = 0.1;
	plane.Ksb = 0.1;
	plane.n = 5.0;
	plane.Iar = 1.0;
	plane.Iag = 1.0;
	plane.Iab = 1.0;
	plane.Isr = 1.0;
	plane.Isg = 1.0;
	plane.Isb = 1.0;

	cube.Kar = 0.2;
	cube.Kag = 0.2;
	cube.Kab = 0.2;
	cube.Kdr = 0.4;
	cube.Kdg = 0.4;
	cube.Kdb = 0.4;
	cube.Ksr = 0.1;
	cube.Ksg = 0.1;
	cube.Ksb = 0.1;
	cube.n = 5.0;
	cube.Iar = 1.0;
	cube.Iag = 1.0;
	cube.Iab = 1.0;
	cube.Isr = 1.0;
	cube.Isg = 1.0;
	cube.Isb = 1.0;

	cone.Kar = 0.3;
	cone.Kag = 0.3;
	cone.Kab = 0.02;
	cone.Kdr = 0.8;
	cone.Kdg = 0.8;
	cone.Kdb = 0.1;
	cone.Ksr = 0.1;
	cone.Ksg = 0.1;
	cone.Ksb = 0.1;
	cone.n = 5.0;
	cone.Iar = 1.0;
	cone.Iag = 1.0;
	cone.Iab = 1.0;
	cone.Isr = 1.0;
	cone.Isg = 1.0;
	cone.Isb = 1.0;

	

	trace();
    	
    glutSwapBuffers();
}
void reshape (int w, int h) {
	/* set the viewport */
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);

	/* Matrix for projection transformation */
	glMatrixMode (GL_PROJECTION); 

	/* replaces the current matrix with the identity matrix */
	glLoadIdentity ();

	/* Define a 2d orthographic projection matrix */
	gluOrtho2D (0.0, (GLdouble) w, 0.0, (GLdouble) h);
}

/*******************************************************************/

int main(int argc, char** argv) {
		/* deal with any GLUT command Line options */
    glutInit(&argc, argv);

	/* create an output window */
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
	glutInitWindowSize(800, 800);

	/* set the name of the window and try to create it */
    glutCreateWindow("CS 488 - Sample");

	/* specify clear values for the color buffers */
	glClearColor (0.0, 0.0, 0.0, 1.0);

    /* Receive keyboard inputs */


    /* assign the display function */
    glutDisplayFunc(display);

	/* assign the idle function */
    glutIdleFunc(display);

    /* sets the reshape callback for the current window */
	glutReshapeFunc(reshape);

    /* enters the GLUT event processing loop */
    glutMainLoop();

    return (EXIT_SUCCESS);
}



