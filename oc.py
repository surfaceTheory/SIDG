#-*- coding: utf-8 -*-
from __future__ import division  
import os, time,sys, util
import coopr.environ
import coopr.pyomo as pym, coopr.opt as opt, numpy as np
#HASM4: two of the main equations of surface theory are used
#HASM5: based on HASM4,add the 3rd equation 
#HASM6: based on HASM5, while boundary is not fixed
mxd=None #None will no export jpg
bndFile=None                                #None will no clip 
symLyr="F:\Py4ArcGIS\OC\gauss.img.lyr"#necessary when above two is provided

def HASM(shpfile, field, cellsize,IniRasFile=None, reTol=1e-5,
          options='HASM6',solver='cplex', out_raster=None):
#get the initial values for a surface#
 if (not options=='HASM6') and (IniRasFile is None):
  print 'initial values is needed for '+ options;   return 0
 import util 
 F0, Zs,RCix,extent, r, c= util.input(shpfile, field,cellsize,IniRasFile=IniRasFile,bndFile=bndFile)
 print solver,',',options, ",Beigin "+options+" modelling:"+str(time.ctime()); t1=time.time()   
 #set up HASM model and declare necessary variables
 md = pym.ConcreteModel('HASMmodel');k=0  #declare a model for HASM application 
 md.F=pym.Var(range(r), range(c)) #declare a two dimension surface  
 
 #grid where there is a sample, is set with the sample's value directly
 #md.constraints = pym.ConstraintList()
 for i in range(len(Zs)):  #Zs[pts*3], rcIx[pts*2]
  md.F[RCix[i,0],RCix[i,1]].value=Zs[i];             md.F[RCix[i,0],RCix[i,1]].fixed=True     #left 
  #md.constraints.add(md.F[RCix[i,0],RCix[i,1]] - Zs[i]==0) 
 #fix the four boundary side
 if  not (options=='HASM6'):  #boundary values is not fixed in HASM6
 #fix the four boundary side
  for i in range(r):
   md.F[i,0].value=F0[i,0];             md.F[i,0].fixed=True     #left 
   md.F[i,c-1].value=F0[i,c-1]; md.F[i,c-1].fixed=True #right 
  for j in range(c):
   md.F[0,j].value=F0[0,j];              md.F[0,j].fixed=True      #bottom 
   md.F[r-1,j].value=F0[r-1,j];  md.F[r-1,j].fixed=True #up  
 rere0=1e38
 print solver, ',', r,  '*' , c, ', reTol=', reTol    
 while(1): #iterate until HASM meet the precision 
  t11=time.time()
  if k>0:  md.del_component('smooth')
  md.smooth=pym.Objective(expr=H456_rule(md,F0,cellsize,options),sense=pym.minimize) 
  #md.smooth.pprint()  
  # Create a solver plugin, other solvers can also be tested here
  optimizer = opt.SolverFactory(solver);   instance = md.create()   
  results = optimizer.solve(instance);     instance.load(results)  
  F1=model2matrix(md,r,c)
  rere=((F1-F0)**2).sum()/(F0**2).sum(); k+=1 #relative residuals  
  print k, ',relative Residual||x(n+1)-x(n)||2/||x(n)||2=,', rere,', ', np.round(time.time()-t11),' seconds'  
  if rere>rere0: 
   print 'Divergence happened, previous results returned.'; 
   print '#'*50  
   F1=F0; break
  else:    rere0=rere  
  
  if rere<=reTol:   break
  else:             F0=F1

 return F1, extent
 
def H456_rule(md,F0,cellsize,options):  
 #minimize the sum of square of variance of all surface grids:
 RHS1,RHS2,RHS3=getRHSofHASM(F0,cellsize, options)      
 #RHS(right hand side) of first and second fundamental equation in HASM
 #minimize the sum of square of the surface theory equations for all grids: 
 r,c=F0.shape; ex=0  
 for i in range(1,r-1):
  for j in range(1,c-1):       
   ex+=(md.F[i+1,j]-2*md.F[i,j]+md.F[i-1,j]-np.float(RHS1[i,j]))**2  #1st equation
   ex+=(md.F[i,j+1]-2*md.F[i,j]+md.F[i,j-1]-np.float(RHS2[i,j]))**2  #2nd equation
   if not (options=='HASM4'):
    ex+=(md.F[i+1,j+1]-md.F[i-1,j+1]+md.F[i-1,j-1]-md.F[i+1,j-1]
         -np.float(RHS3[i,j]))**2  #3rd equation  
 del RHS1,RHS2,RHS3                                         
 return ex  

 
def getRHSofHASM(Ma,h, options):
#  RHS1=np.zeros_like(Ma)#RHS of first fundamental equation in HASM 
#  RHS2=np.zeros_like(Ma)#RHS of second fundamental equation in HASM
#  RHS3=np.zeros_like(Ma)#RHS of third fundamental equation in HASM
 r,c=Ma.shape 
 fx,fy=np.gradient(Ma,h);fxx,fxy=np.gradient(fx,h);fxy,fyy=np.gradient(fy,h)
 
 E=1+fx**2;  F=fx*fy; G=1+fy**2; eg=2*(E + G - 1);  eg2=(eg/2)**0.5
 L=fxx/eg2;  N=fyy/eg2
 
 Ex,Ey=np.gradient(E,h);  Fx, Fy=np.gradient(F,h); Gx,Gy=np.gradient(G,h) 
 gm111=(G*Ex - 2*F*Fx + F*Ey)/eg;   gm221=(2*G*Fy - G*Gx - F*Gy)/eg
 gm112=(2*E*Fx - E*Ey - F*Ex)/eg;   gm222=(E*Gy - 2*F*Fy + F*Gx)/eg
 
 if not (options == 'HASM4'):
  M=fxy/eg2;  gama121=(G*Ey-F*Gx)/eg; gama122=(E*Gx-F*Ey)/eg   
 
 for i in range(1,r-1):
  for j in range(1, c-1):
   fxx[i,j]=(0.5*gm111[i,j]*(Ma[i+1,j]-Ma[i-1,j])/h
            + 0.5*gm112[i,j]*(Ma[i,j+1]-Ma[i,j-1])/h+ L[i,j]/eg2[i,j])*h**2
   fyy[i,j]=(0.5*gm221[i,j]*(Ma[i+1,j]-Ma[i-1,j])/h
            + 0.5*gm222[i,j]*(Ma[i,j+1]-Ma[i,j-1])/h+ N[i,j]/eg2[i,j])*h**2   
   if not (options=='HASM4'):             
    fxy[i,j]=(0.5*gama121[i,j]*(Ma[i+1,j]-Ma[i-1,j])/h
    +0.5*gama122[i,j]*(Ma[i,j+1]-Ma[i,j-1])/h+M[i,j]/eg2[i,j])*4*h**2    
 del E,F,G,L,N,fx,fy, Ex,Ey,Fx,Fy,Gx,Gy, eg,eg2
 del gm111,gm112,gm221,gm222
 if not (options=='HASM4'):
  del M,gama121, gama122 
 return fxx,fyy,fxy 


def model2matrix22(md,r,c):
 Fvalues=[md.F[i,j].value for i in range(r) for j in range(c)]
 W=np.fromiter(Fvalues, np.float);    W.shape=(r,c);    return W

def model2matrix(md,r,c):
    Fvalues2=[md.F[i,j].value for i in range(r) for j in range(c)]
    Fvalues=[item[2] for item in Fvalues2 ]
    W=np.fromiter(Fvalues, np.float);   
    W.shape=(r,c);    
    return W 
 
if __name__=='__main__': 
    data='gauss5uf5';   solver='knitro'
    shpfile='c:\\HASMCode\\gauss\\'+data+'.shp';    field='GRID_CODE'; cellsize=6.0/30 
    HASM(shpfile=shpfile, field=field, cellsize=cellsize,IniRasFile=None, reTol=1e-5,
          options='HASM6',solver=solver, out_raster='c:\\hasmcode\\gass.img')
     
