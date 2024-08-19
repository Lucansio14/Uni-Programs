# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:47:43 2019
Surface Area and Volume Calculator for 3D Geometric Figures
@author: Lucas Romero Fernández
"""
#All units in International System units (S.I.) (for example) to avoid confusion.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
#main_program
FigAvailable={1:"Parallelepiped",2:"Sphere",3:"Circular cylinder",4:"Arbitrary cylinder (no figure due to arbitrariness)",5:"Right circular cone",6:"Pyramid of arbitrary base (only volume and no figure due to arbitrariness)",7:"Spherical cap (in a sphere)",8:"Frustrum of a right circular cone",9:"Spherical triangle (no figure)",10:"Ring torus",11:"Ellipsoid (only volume)",12:"Paraboloid of revolution (only volume)"}
print("*Beware that all units of the variables have to be cohesive (for example, all expressed in the International System units (S.I.)) to avoid confusion and all values given by the user will be considered positive (except angles), all of this with the purpose of obtaining coherent results. Also, for ease of calculus, a vertex of the figure (or, in specific cases, the center of the figure) always coincides with the origin of coordinates (0,0,0).*")
print("")
print("3D geometric figures available: ")
print("")
print(FigAvailable)
print("")
while True:
    FigSelector=abs(int(input("Type the number/integer corresponding to the position of the desired 3D geometric figure of the list above: ")))
    if FigSelector>12:
        print("")
        print("Type an integer that represents the position of the desired 3D geometric figure of the list above, please.")
        print("")
        continue
    print("")
    print("Figure selected:",FigAvailable[FigSelector])
    print("")
    #Surface area and volume of a parallelepiped
    if FigSelector==1:
        side1Para=abs(float(input("Length of side 1 of the parallelepiped? ")))
        print("")
        side2Para=abs(float(input("Length of side 2 of the parallelepiped? ")))
        print("")
        side3Para=abs(float(input("Length of side 3 of the parallelepiped? ")))
        print("")
        angleParaDeg1=float(input("Angle (in degrees) between side 1 and side 2 of the parallelepiped? "))
        print("")
        angleParaDeg2=float(input("Angle (in degrees) between side 2 and side 3 of the parallelepiped? "))
        print("")
        angleParaRad1=np.radians(angleParaDeg1)
        angleParaRad2=np.radians(angleParaDeg2)
        surfareaPara=2*(side1Para*side2Para*np.sin(angleParaRad1)+side2Para*side3Para*np.sin(angleParaRad2)+side1Para*side3Para)
        volPara=side1Para*side2Para*side3Para*np.sin(angleParaRad2)
        pointsPara1=np.array([[0,0,0],[side1Para,0,0],[side1Para+side2Para*np.cos(angleParaRad1),side2Para*np.sin(angleParaRad1),0],[side2Para*np.cos(angleParaRad1),side2Para*np.sin(angleParaRad1),0],[0,side3Para*np.cos(angleParaRad2),side3Para*np.sin(angleParaRad2)],[side1Para,side3Para*np.cos(angleParaRad2),side3Para*np.sin(angleParaRad2)],[side1Para+side2Para*np.cos(angleParaRad1),side2Para*np.sin(angleParaRad1)+side3Para*np.cos(angleParaRad2),side3Para*np.sin((angleParaRad2))],[side2Para*np.cos(angleParaRad1),side2Para*np.sin(angleParaRad1)+side3Para*np.cos(angleParaRad2),side3Para*np.sin(angleParaRad2)]])
        print("The parallelepiped of sides of length ",side1Para,", ",side2Para," and ",side3Para," and angle between the sides 1 and 2 (in degrees) of ",angleParaDeg1,"and angle between the sides 2 and 3 (in degrees) of ",angleParaDeg2," has a surface area of value ",surfareaPara," and a volume of value ",volPara,".",sep="")
        P=[[2.06498904e-01,-6.30755443e-07,1.07477548e-03],[1.61535574e-06,1.18897198e-01,7.85307721e-06],[7.08353661e-02,4.48415767e-06,2.05395893e-01]]
        Z=np.zeros((8,3))
        for i in range(8):
            Z[i,:]=np.dot(pointsPara1[i,:],P)
        Z=10.0*Z
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        pointsPara2=np.array([[Z[0]],[Z[1]],[Z[2]],[Z[3]],[Z[4]],[Z[5]],[Z[6]],[Z[7]]])
        ax.scatter3D(Z[:,0],Z[:,1],Z[:,2])
        vertices=[[Z[0],Z[1],Z[2],Z[3]],[Z[4],Z[5],Z[6],Z[7]],[Z[0],Z[1],Z[5],Z[4]],[Z[2],Z[3],Z[7],Z[6]],[Z[1],Z[2],Z[6],Z[5]],[Z[4],Z[7],Z[3],Z[0]]]
        ax.add_collection3d(Poly3DCollection(vertices,facecolors="black",linewidths=1,edgecolors="white",alpha=0.4))
        ax.set_aspect("equal","box")
        plt.title("Parallelepiped")
        plt.grid()
        plt.tight_layout()
        plt.show()
        print("Vertex coordinates: ",pointsPara2,sep="")
        break
    #Surface area and volume of a sphere
    if FigSelector==2:
        radiSphere=abs(float(input("Radius of the sphere? ")))
        print("")
        surfareaSphere=4*np.pi*(radiSphere)**2
        volSphere=(4/3)*np.pi*(radiSphere)**3
        print("The sphere of radius ",radiSphere," has a surface area of value ",surfareaSphere," and a volume of value ",volSphere,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        phi=np.linspace(0,2*np.pi,100)
        theta=np.linspace(0,np.pi,100)
        theta,phi=np.meshgrid(theta,phi)
        x=radiSphere*np.cos(phi)*np.sin(theta)
        y=radiSphere*np.sin(phi)*np.sin(theta)
        z=radiSphere*np.cos(theta)
        ax.plot_wireframe(x,y,z,rstride=5,cstride=5,color="black")
        plt.plot(0,0,0,marker="o",color="red")
        ax.set_aspect("equal","box")
        plt.title("Sphere")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Surface area and volume of a circular cylinder
    if FigSelector==3:
        radiCircCyl=abs(float(input("Radius of the circular cylinder? ")))
        print("")
        inclheightCircCyl=abs(float(input("Inclined height (or the height of the cylinder if it is right) of the circular cylinder? ")))
        print("")
        angleCircCylDeg=float(input("Angle (in degrees) between the base and the inclined height of the circular cylinder? "))
        print("")
        angleCircCylRad=np.radians(angleCircCylDeg)
        surfareaCircCyl=2*np.pi*radiCircCyl*inclheightCircCyl+2*np.pi*(radiCircCyl)**2
        volCircCyl=np.pi*inclheightCircCyl*np.sin(angleCircCylRad)*(radiCircCyl)**2
        print("The circular cylinder of radius ",radiCircCyl,", inclined height (or just the height if it is a right circular cylinder) of ",inclheightCircCyl," and angle between the base and the inclined height (in degrees) of ",angleCircCylDeg," has a surface area of value ",surfareaCircCyl," and a volume of value ",volCircCyl,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        theta=np.linspace(0,2*np.pi,100)
        u=np.linspace(0,inclheightCircCyl*np.sin(angleCircCylRad),100)
        u,theta=np.meshgrid(u,theta)
        x=radiCircCyl*np.cos(theta)+u*np.cos(angleCircCylRad)
        y=radiCircCyl*np.sin(theta)
        z=u*np.sin(angleCircCylRad)
        ax.plot_wireframe(x,y,z,rstride=5,cstride=5,color="black")
        plt.plot(0,0,0,marker="o",color="red")
        plt.plot(inclheightCircCyl*np.sin(angleCircCylRad)*np.sin(np.pi/2-angleCircCylRad),0,inclheightCircCyl*np.sin(angleCircCylRad)*np.cos(np.pi/2-angleCircCylRad),marker="o",color="red")
        ax.set_aspect("equal","box")
        plt.title("Circular cylinder")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Surface area and volume of an arbitrary cylinder
    if FigSelector==4:
        crosssectareaArbCyl=abs(float(input("Cross-sectional area of the arbitrary cylinder? ")))
        print("")
        periArbCyl=abs(float(input("Perimeter of the cover and base of the arbitrary cylinder? ")))
        print("")
        inclheightArbCyl=abs(float(input("Inclined height of the arbitrary cylinder? ")))
        print("")
        angleArbCylDeg=float(input("Angle (in degrees) between the base and the inclined height of the arbitrary cylinder? "))
        print("")
        angleArbCylRad=np.radians(angleArbCylDeg)
        surfareaArbCyl=periArbCyl*inclheightArbCyl*np.sin(angleArbCylRad)+2*crosssectareaArbCyl
        volArbCyl=crosssectareaArbCyl*inclheightArbCyl*np.sin(angleArbCylRad)
        print("The arbitrary cylinder of cross-sectional area of ",crosssectareaArbCyl,", inclined height of ",inclheightArbCyl,", perimeter of the cover and base of ",periArbCyl," and angle between the base and the inclined height (in degrees) of ",angleArbCylDeg," has a surface area of value ",surfareaArbCyl," and a volume of value ",volArbCyl,".",sep="")
        break
    #Surface area and volume of a right circular cone
    if FigSelector==5:
        radiCircCone=abs(float(input("Radius of the right circular cone? ")))
        print("")
        heightCircCone=float(input("Height of the right circular cone? "))
        print("")
        surfareaCircCone=np.pi*radiCircCone*np.sqrt(radiCircCone**2+heightCircCone**2)+2*np.pi*(radiCircCone)**2
        volCircCone=(1/3)*np.pi*abs(heightCircCone)*(radiCircCone)**2
        print("The right circular cone of radius ",radiCircCone," and height of ",heightCircCone," has a surface area of value ",surfareaCircCone," and a volume of value ",volCircCone,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        theta=np.linspace(0,2*np.pi,100)
        u=np.linspace(0,heightCircCone,100)
        u,theta=np.meshgrid(u,theta)
        x=(radiCircCone-(radiCircCone/heightCircCone)*u)*np.cos(theta)
        y=(radiCircCone-(radiCircCone/heightCircCone)*u)*np.sin(theta)
        z=u
        ax.plot_wireframe(x,y,z,rstride=5,cstride=5,color="black")
        plt.plot(0,0,0,marker="o",color="red")
        ax.set_aspect("equal","box")
        plt.title("Right circular cone")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Volume of a pyramid of arbitrary base
    if FigSelector==6:
        areaPyr=abs(float(input("Area of the base of the pyramid of arbitrary base? ")))
        print("")
        heightPyr=abs(float(input("Height of the pyramid of arbitrary base? ")))
        print("")
        volPyr=(1/3)*areaPyr*heightPyr
        print("The pyramid of arbitrary base of area of the base of ",areaPyr," and height of ",heightPyr," has a volume of value ",volPyr,".",sep="")
        break
    #Surface area and volume of a spherical cap (in a sphere)
    if FigSelector==7:
        radiSphCap=abs(float(input("Radius of the sphere that includes the spherical cap? ")))
        print("")
        heightSphCap=abs(float(input("Height of the spherical cap? ")))
        if heightSphCap>2*radiSphCap:
            sys.exit("The height of the spherical cap cannot be bigger than two times the radius (the diameter, in other words) of the sphere that includes the spherical cap...")
        print("")
        surfareaSphCap=2*np.pi*radiSphCap*heightSphCap
        volSphCap=(1/3)*np.pi*(heightSphCap**2)*(3*radiSphCap-heightSphCap)
        print("The spherical cap of radius of the sphere that includes the spherical cap of ",radiSphCap," and height of the spherical cap of ",heightSphCap," has a surface area of value ",surfareaSphCap," and a volume of value ",volSphCap,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        phi=np.linspace(0,2*np.pi,100)
        theta=np.linspace(0,np.arccos(1-(heightSphCap/radiSphCap)),100)
        theta,phi=np.meshgrid(theta,phi)
        x=radiSphCap*np.cos(phi)*np.sin(theta)
        y=radiSphCap*np.sin(phi)*np.sin(theta)
        z=radiSphCap*np.cos(theta)
        ax.plot_surface(x,y,z,rstride=5,cstride=5,color="black",edgecolors="white")
        ax.set_aspect("equal","box")
        plt.title("Spherical cap")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Surface area and volume of a frustrum of a right circular cone
    if FigSelector==8:
        radiFrustRightCircCone1=abs(float(input("Radius of the smaller circle of the frustrum of a right circular cone? ")))
        print("")
        radiFrustRightCircCone2=abs(float(input("Radius of the bigger circle of the frustrum of a right circular cone? ")))
        print("")
        heightFrustRightCircCone=float(input("Height of the frustrum of a right circular cone? "))
        print("")
        surfareaFrustRightCircCone=np.pi*(radiFrustRightCircCone1+radiFrustRightCircCone2)*np.sqrt((heightFrustRightCircCone)**2+(radiFrustRightCircCone1-radiFrustRightCircCone2)**2)+np.pi*(radiFrustRightCircCone1)**2+np.pi*(radiFrustRightCircCone2)**2
        volFrustRightCircCone=(1/3)*np.pi*abs(heightFrustRightCircCone)*(radiFrustRightCircCone2**2+radiFrustRightCircCone1*radiFrustRightCircCone2+radiFrustRightCircCone1**2)
        print("The frustrum of a right circular cone of radius of the smaller circle of ",radiFrustRightCircCone1,", radius of the bigger circle of ",radiFrustRightCircCone2," and height of ",heightFrustRightCircCone," has a surface area of value ",surfareaFrustRightCircCone," and a volume of value ",volFrustRightCircCone,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        theta=np.linspace(0,2*np.pi,100)
        u=np.linspace(0,heightFrustRightCircCone,100)
        u,theta=np.meshgrid(u,theta)
        x=(radiFrustRightCircCone2-(radiFrustRightCircCone1/heightFrustRightCircCone)*u)*np.cos(theta)
        y=(radiFrustRightCircCone2-(radiFrustRightCircCone1/heightFrustRightCircCone)*u)*np.sin(theta)
        z=u
        ax.plot_wireframe(x,y,z,rstride=5,cstride=5,color="black")
        plt.plot(0,0,0,marker="o",color="red")
        plt.plot(0,0,heightFrustRightCircCone,marker="o",color="red")
        ax.set_aspect("equal","box")
        plt.title("Frustrum of a right circular cone")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Area of a spherical triangle
    if FigSelector==9:
        radiSphTri=abs(float(input("Radius of the sphere where the spherical triangle is on? ")))
        print("")
        angleSphTriDeg1=float(input("Angle 1 (in degrees) of the spherical triangle? "))
        print("")
        angleSphTriDeg2=float(input("Angle 2 (in degrees) of the spherical triangle? "))
        print("")
        angleSphTriDeg3=float(input("Angle 3 (in degrees) of the spherical triangle? "))
        print("")
        angleSphTriRad1=np.radians(angleSphTriDeg1)
        angleSphTriRad2=np.radians(angleSphTriDeg2)
        angleSphTriRad3=np.radians(angleSphTriDeg3)
        areaSphTri=abs((angleSphTriRad1+angleSphTriRad2+angleSphTriRad3-np.pi)*(radiSphTri)**2)
        print("The spherical triangle of angle 1 (in degrees) of ",angleSphTriDeg1,", angle 2 (in degrees) of ",angleSphTriDeg2,", angle 3 (in degrees) of ",angleSphTriDeg3," and radius of the sphere where the spherical triangle is on of ",radiSphTri," has an area of value ",areaSphTri,".",sep="")
        break
    #Surface area and volume of a ring torus
    if FigSelector==10:
        radiintTorus=abs(float(input("Interior radius of the ring torus? ")))
        print("")
        radiextTorus=abs(float(input("Exterior radius of the ring torus? ")))
        print("")
        if radiintTorus>radiextTorus:
            print("It looks like you have swapped the values of the two radii, let me fix it...")
            print("")
            change=radiintTorus
            radiintTorus=radiextTorus
            radiextTorus=change
        surfareaTorus=(np.pi**2)*(radiextTorus**2-radiintTorus**2)
        volTorus=(1/4)*(np.pi**2)*(radiintTorus+radiextTorus)*((radiextTorus-radiintTorus)**2)
        print("The ring torus of interior radius of ",radiintTorus," and exterior radius of ",radiextTorus," has a surface area of value ",surfareaTorus," and a volume of value ",volTorus,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        phi=np.linspace(0,2*np.pi,100)
        theta=np.linspace(0,2*np.pi,100)
        theta,phi=np.meshgrid(theta,phi)
        r=(radiextTorus-radiintTorus)/2
        R=radiintTorus+r
        x=(R+r*np.cos(theta))*np.cos(phi)
        y=(R+r*np.cos(theta))*np.sin(phi)
        z=r*np.sin(theta)
        ax.plot_surface(x,y,z,rstride=5,cstride=5,color="black",edgecolors="white")
        ax.view_init(36,26)
        plt.plot(0,0,0,marker="o",color="red")
        ax.set_aspect("equal","box")
        plt.title("Ring torus")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Volume of an ellipsoid
    if FigSelector==11:
        semiaxisxEll=abs(float(input("Length of the semi-axis x of the ellipsoid? ")))
        print("")
        semiaxisyEll=abs(float(input("Length of the semi-axis y of the ellipsoid? ")))
        print("")
        semiaxiszEll=abs(float(input("Length of the semi-axis z of the ellipsoid? ")))
        print("")
        volEll=(4/3)*np.pi*semiaxisxEll*semiaxisyEll*semiaxiszEll
        print("The ellipsoid of semi-axis x of length ",semiaxisxEll,", semi-axis y of length ",semiaxisyEll," and semi-axis z of length ",semiaxiszEll," has a volume of value ",volEll,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        phi=np.linspace(0,2*np.pi,100)
        theta=np.linspace(0,np.pi,100)
        theta,phi=np.meshgrid(theta,phi)
        x=semiaxisxEll*np.cos(phi)*np.sin(theta)
        y=semiaxisyEll*np.sin(phi)*np.sin(theta)
        z=semiaxiszEll*np.cos(theta)
        ax.plot_wireframe(x,y,z,rstride=5,cstride=5,color="black")
        plt.plot(0,0,0,marker="o",color="red")
        ax.set_aspect("equal","box")
        plt.title("Ellipsoid")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Volume of a paraboloid of revolution
    if FigSelector==12:
        maxminheightParRev=float(input("Maximum (or minimum) height of the paraboloid of revolution? "))
        print("")
        radibaseParRev=abs(float(input("Radius of the base circle of the paraboloid of revolution? ")))
        print("")
        volParRev=(1/2)*np.pi*abs(maxminheightParRev)*(radibaseParRev)**2
        print("The paraboloid of revolution of maximum height of ",maxminheightParRev," and radius of the base circle of ",radibaseParRev," has a volume of value ",volParRev,".",sep="")
        fig=plt.figure()
        ax=fig.add_subplot(projection="3d")
        phi=np.linspace(0,2*np.pi,100)
        u=np.linspace(0,maxminheightParRev,100)
        u,phi=np.meshgrid(u,phi)
        x=radibaseParRev*np.sqrt(u/maxminheightParRev)*np.cos(phi)
        y=radibaseParRev*np.sqrt(u/maxminheightParRev)*np.sin(phi)
        z=u
        ax.plot_surface(x,y,z,rstride=5,cstride=5,color="black",edgecolors="white")
        ax.set_aspect("equal","box")
        plt.title("Paraboloid of revolution")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
