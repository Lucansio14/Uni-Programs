# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:47:43 2019
Perimeter and Area (and miscellaneous) Calculator for Flat Geometric Figures
@author: Lucas Romero Fernández
"""
#All units in International System units (S.I.) (for example) to avoid confusion.
#main_program
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib.patches as pltpat
import sympy as sy
from sympy.solvers import solve
class Shape:
    def __init__(self, color='Black'):
        self.color=color
    def get_color(self):
        return self.color
    def draw(self):
        pass
class Parallelogram(Shape):
    def __init__(self,points):
        super().__init__()
        self.points=points
    def draw(self):
        parallelogram=pltpat.Polygon(self.points,color=self.color,fill=False)
        fig,ax=plt.subplots()
        x,y=self.points.T
        plt.scatter(x,y,marker="o",color=self.color)
        ax.add_patch(parallelogram)
        ax.set_aspect("equal","box")
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Triangle(Shape):
    def __init__(self,points):
        super().__init__()
        self.points=points
    def draw(self):
        triangle=pltpat.Polygon(self.points,color=self.color,fill=False)
        fig,ax=plt.subplots()
        x,y=self.points.T
        plt.scatter(x,y,marker="o",color=self.color)
        ax.add_patch(triangle)
        ax.set_aspect("equal","box")
        plt.title(type(self).__name__)  
        plt.grid()
        plt.tight_layout()
        plt.show()
class Trapezoid(Shape):
    def __init__(self,points):
        super().__init__()
        self.points=points
    def draw(self):
        trapezoid=pltpat.Polygon(self.points,color=self.color,fill=False)
        fig,ax=plt.subplots()
        x,y=self.points.T
        plt.scatter(x,y,marker="o",color=self.color)
        ax.add_patch(trapezoid)
        ax.set_aspect("equal","box")
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Regular_Polygon(Shape):
    def __init__(self,radius,numVertices,circu,radius2=0):
        super().__init__()
        self.radius=radius
        self.numVertices=numVertices
        self.circu=circu
        self.radius2=radius2
    def draw(self):
        regpolygon=pltpat.RegularPolygon((0,0),self.numVertices,radius=self.radius,color=self.color,fill=False)
        if self.circu==0:
            circle=plt.Circle((0,0),self.radius,color=self.color,fill=False,alpha=0.7)
        if self.circu==1:
            circle=plt.Circle((0,0),self.radius2,color=self.color,fill=False,alpha=0.7)
        fig,ax=plt.subplots()
        plt.plot(0,0,marker="o",color=self.color)
        if self.circu==0:
            plt.axhline(0,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axhline(-self.radius,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axhline(self.radius,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axvline(0,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axvline(-self.radius,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axvline(self.radius,color=self.color,linestyle="dashed",alpha=0.4)
        if self.circu==1:
            plt.axhline(0,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axhline(-self.radius2,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axhline(self.radius2,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axvline(0,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axvline(-self.radius2,color=self.color,linestyle="dashed",alpha=0.4)
            plt.axvline(self.radius2,color=self.color,linestyle="dashed",alpha=0.4)
        ax.add_patch(regpolygon)
        ax.add_patch(circle)
        ax.set_aspect("equal","box")
        if self.circu==0:
            plt.xlim(-self.radius-1,self.radius+1)
            plt.ylim(-self.radius-1,self.radius+1)
        if self.circu==1:
            plt.xlim(-self.radius2-1,self.radius2+1)
            plt.ylim(-self.radius2-1,self.radius2+1)
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Circle(Shape):
    def __init__(self,radius):
        super().__init__()
        self.radius=radius
    def draw(self):
        circle=plt.Circle((0,0),self.radius,color=self.color,fill=False)
        fig,ax=plt.subplots()
        plt.plot(0,0,marker="o",color=self.color)
        plt.axhline(0,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axhline(-self.radius,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axhline(self.radius,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(0,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(-self.radius,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(self.radius,color=self.color,linestyle="dashed",alpha=0.4)
        ax.add_patch(circle)
        ax.set_aspect("equal","box")
        plt.xlim(-self.radius-1,self.radius+1)
        plt.ylim(-self.radius-1,self.radius+1)
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Circular_Sector(Shape):
    def __init__(self,radius,theta):
        super().__init__()
        self.radius=radius
        self.theta=theta
    def draw(self):
        circularsector=pltpat.Wedge((0,0),self.radius,0,self.theta,color=self.color,fill=False)
        circle=plt.Circle((0,0),self.radius,color=self.color,fill=False,alpha=0.5)
        fig,ax=plt.subplots()
        plt.plot(0,0,marker="o",color=self.color)
        plt.axhline(0,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axhline(-self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axhline(self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axvline(0,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axvline(-self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axvline(self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        ax.add_patch(circularsector)
        ax.add_patch(circle)
        ax.set_aspect("equal","box")
        plt.xlim(-self.radius-1,self.radius+1)
        plt.ylim(-self.radius-1,self.radius+1)
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Circle_Inscribed_In_Triangle(Shape):
    def __init__(self,points,radius,incenterx,incentery):
        super().__init__()
        self.points=points
        self.radius=radius
        self.incenterx=incenterx
        self.incentery=incentery
    def draw(self):
        triangle=pltpat.Polygon(self.points,color=self.color,fill=False)
        circle=plt.Circle((self.incenterx,self.incentery),self.radius,color=self.color,fill=False)
        fig,ax=plt.subplots()
        x,y=self.points.T
        plt.scatter(x,y,marker="o",color=self.color)
        plt.plot(self.incenterx,self.incentery,marker="o",color=self.color)
        plt.axhline(self.incentery,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(self.incenterx,color=self.color,linestyle="dashed",alpha=0.4)
        ax.add_patch(triangle)
        ax.add_patch(circle)
        ax.set_aspect("equal","box")
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Circle_Circumscribed_By_Triangle(Shape):
    def __init__(self,points,radius,circumcenterx,circumcentery):
        super().__init__()
        self.points=points
        self.radius=radius
        self.circumcenterx=circumcenterx
        self.circumcentery=circumcentery
    def draw(self):
        triangle=pltpat.Polygon(self.points,color=self.color,fill=False)
        circle=plt.Circle((self.circumcenterx,self.circumcentery),self.radius,color=self.color,fill=False)
        fig,ax=plt.subplots()
        x,y=self.points.T
        plt.scatter(x,y,marker="o",color=self.color)
        plt.plot(self.circumcenterx,self.circumcentery,marker="o",color=self.color)
        plt.axhline(self.circumcentery,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(self.circumcenterx,color=self.color,linestyle="dashed",alpha=0.4)
        ax.add_patch(triangle)
        ax.add_patch(circle)
        ax.set_aspect("equal","box")
        plt.xlim(self.circumcenterx-self.radius-1,self.circumcenterx+self.radius+1)
        plt.ylim(self.circumcentery-self.radius-1,self.circumcentery+self.radius+1)
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Circle_Segment(Shape):
    def __init__(self,points,radius,theta):
        super().__init__()
        self.points=points
        self.radius=radius
        self.theta=theta
    def draw(self):
        circularsector=pltpat.Wedge((0,0),self.radius,(180-self.theta)/2,((180-self.theta)/2)+self.theta,color=self.color)
        circle=plt.Circle((0,0),self.radius,color=self.color,fill=False,alpha=0.5)
        triangle=pltpat.Polygon(self.points,color="white")
        fig,ax=plt.subplots()
        plt.plot(0,0,marker="o",color=self.color)
        plt.axhline(0,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axhline(-self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axhline(self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axvline(0,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axvline(-self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        plt.axvline(self.radius,color=self.color,linestyle="dashed",alpha=0.2)
        ax.add_patch(circularsector)
        ax.add_patch(circle)
        ax.add_patch(triangle)
        ax.set_aspect("equal","box")
        plt.xlim(-self.radius-1,self.radius+1)
        plt.ylim(-self.radius-1,self.radius+1)
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
class Ellipse(Shape):
    def __init__(self,width,height):
        super().__init__()
        self.width=width
        self.height=height
    def draw(self):
        ellipse=pltpat.Ellipse((0,0),self.width,self.height,color=self.color,fill=False)
        fig,ax=plt.subplots()
        plt.plot(0,0,marker="o",color=self.color)
        plt.axhline(0,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axhline(-(1/2)*self.height,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axhline((1/2)*self.height,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(0,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline(-(1/2)*self.width,color=self.color,linestyle="dashed",alpha=0.4)
        plt.axvline((1/2)*self.width,color=self.color,linestyle="dashed",alpha=0.4)
        ax.add_patch(ellipse)
        ax.set_aspect("equal","box")
        plt.xlim(-(1/2)*self.width-1,(1/2)*self.width+1)
        plt.ylim(-(1/2)*self.height-1,(1/2)*self.height+1)
        plt.title(type(self).__name__)
        plt.grid()
        plt.tight_layout()
        plt.show()
FigAvailable={1:"Parallelogram",2:"Triangle",3:"Trapezoid",4:"Regular polygon",5:"Circle",6:"Circular sector",7:"Circle inscribed in a triangle (to calculate radius of that circle)",8:"Circle circumscribed by a triangle (to calculate radius of that circle)",9:"Regular polygon inscribed in a circle",10:"Regular polygon circumscribed by a circle",11:"Circle segment (specifically, the region of a disk which is *cut off* from the rest of the disk by a straight line)",12:"Parabolic segment",13:"Ellipse"}
print("*Beware that all units of the variables have to be cohesive (for example, all expressed in the International System units (S.I.)) to avoid confusion and all values given by the user will be considered positive (except angles), all of this with the purpose of obtaining coherent results.*")
print("")
print("Flat geometric figures available: ")
print("")
print(FigAvailable)
print("")
while True:
    FigSelector=abs(int(input("Type the number/integer corresponding to the position of the desired flat geometric figure of the list above: ")))
    if FigSelector>13:
        print("")
        print("Type an integer that represents the position of the desired 3D geometric figure of the list above, please.")
        print("")
        continue
    print("")
    print("Figure selected:",FigAvailable[FigSelector])
    print("")
    #Area and perimeter of a parallelogram
    if FigSelector==1:
        side1Para=abs(float(input("Length of side 1 of the parallelogram? ")))
        print("")
        side2Para=abs(float(input("Length of side 2 of the parallelogram? ")))
        print("")
        angleParaDeg=float(input("Angle (in degrees) between side 1 and side 2 of the parallelogram? "))
        print("")
        angleParaRad=np.radians(angleParaDeg)
        areaPara=side1Para*side2Para*np.sin(angleParaRad)
        periPara=2*side1Para+2*side2Para
        pointsPara=np.array([[0,0],[side1Para*np.cos(angleParaRad),side1Para*np.sin(angleParaRad)],[side1Para*np.cos(angleParaRad)+side2Para,side1Para*np.sin(angleParaRad)],[side2Para,0]])
        print("The parallelogram of sides of length ",side1Para," and ",side2Para," and angle between the two (in degrees) of ",angleParaDeg," has an area of value ",areaPara," and a perimeter of value ",periPara,".",sep="")
        parallelogram=Parallelogram(points=pointsPara)
        parallelogram.color="Black"
        parallelogram.draw()
        print("Vertex coordinates: ",pointsPara,sep="")
        break
    #Area and perimeter of a triangle
    if FigSelector==2:
        side1Tri=abs(float(input("Length of side 1 of the triangle? ")))
        print("")
        side2Tri=abs(float(input("Length of side 2 of the triangle? ")))
        print("")
        side3Tri=abs(float(input("Length of side 3 of the triangle? ")))
        print("")
        angleTriDeg=float(input("Angle (in degrees) between side 1 and side 2 of the triangle? "))
        print("")
        angleTriRad=np.radians(angleTriDeg)
        areaTri=0.5*side1Tri*side2Tri*np.sin(angleTriRad)
        periTri=side1Tri+side2Tri+side3Tri
        x=sy.Symbol("x")
        y=sy.Symbol("y")
        eq1=x**2+y**2-side2Tri**2
        eq2=(x-side1Tri*np.cos(angleTriRad))**2+(y-side1Tri*np.sin(angleTriRad))**2-side3Tri**2
        point3=solve([eq1,eq2],[x,y],dict=True)
        print(point3)
        values=list(point3[1].values())
        pointsTri=np.array([[0,0],[side1Tri*np.cos(angleTriRad),side1Tri*np.sin(angleTriRad)],[values[0],values[1]]])
        print("The triangle of sides of length ",side1Tri,", ",side2Tri," and ",side3Tri," and angle between the first two sides (in degrees) of ",angleTriDeg," has an area of value ",areaTri," and a perimeter of value ",periTri,".",sep="")
        triangle=Triangle(points=pointsTri)
        triangle.color="Black"
        triangle.draw()
        print("Vertex coordinates: ",pointsTri,sep="")
        break
    #Area and perimeter of a trapezoid
    if FigSelector==3:
        base1Tra=abs(float(input("Length of the inferior base of the trapezoid? ")))
        print("")
        heightTra=abs(float(input("Height of the trapezoid? ")))
        print("")
        angleTraDeg1=float(input("Angle 1 (in degrees) of the trapezoid? "))
        print("")
        angleTraDeg2=float(input("Angle 2 (in degrees) of the trapezoid? "))
        print("")
        angleTraRad1=np.radians(angleTraDeg1)
        angleTraRad2=np.radians(angleTraDeg2)
        base2Tra=base1Tra-heightTra*(1/np.tan(angleTraRad1))-heightTra*(1/np.tan(angleTraRad2))
        areaTra=0.5*heightTra*(base1Tra+base2Tra)
        periTra=base1Tra+base2Tra*heightTra*((1/(np.sin(angleTraRad1)))+(1/(np.sin(angleTraRad2))))
        pointsTra=np.array([[0,0],[heightTra*(1/np.tan(angleTraRad1)),heightTra],[base1Tra-heightTra*(1/np.tan(angleTraRad2)),heightTra],[base1Tra,0]])
        print("The trapezoid of sides of length ",base1Tra," and ",base2Tra," and of height ",heightTra," and angles (in degrees) of ",angleTraDeg1," and ",angleTraDeg2," has an area of value ",areaTra," and a perimeter of value ",periTra,".",sep="")
        trapezoid=Trapezoid(points=pointsTra)
        trapezoid.color="Black"
        trapezoid.draw()
        print("Vertex coordinates: ",pointsTra,sep="")
        break
    #Area and perimeter of a regular polygon
    if FigSelector==4:
        NSidesPol=abs(int(input("Number of sides of the regular polygon? ")))
        print("")
        sidePol=abs(float(input("Length of any side of the regular polygon? ")))
        print("")
        areaPol=(1/4)*float(NSidesPol)*(sidePol)**2*((np.cos((np.pi/(float(NSidesPol))))/(np.sin((np.pi)/(float(NSidesPol))))))
        periPol=float(NSidesPol)*sidePol
        radiPol=0.5*sidePol*float(mp.csc(np.pi/NSidesPol))
        print("The regular polygon of ",NSidesPol," sides of length ",sidePol," has an area of value ",areaPol," and a perimeter of value ",periPol,".",sep="")
        regpolygon=Regular_Polygon(radius=radiPol,numVertices=NSidesPol,circu=0)
        regpolygon.color="Black"
        regpolygon.draw()
        break
    #Area and perimeter of a circle
    if FigSelector==5:
        radiCirc=abs(float(input("Radius of the circle? ")))
        print("")
        areaCirc=np.pi*radiCirc**2
        periCirc=2*np.pi*radiCirc
        print("The circle of radius ",radiCirc," has an area of value ",areaCirc," and a perimeter of value ",periCirc,".",sep="")
        circle=Circle(radius=radiCirc)
        circle.color="Black"
        circle.draw()
        break
    #Area and arc length of a circular sector
    if FigSelector==6:
        radiCircSec=abs(float(input("Radius of the circular sector? ")))
        print("")
        angleCircSecDeg=float(input("Angle (in degrees) formed between the two line segments joining the center to the end-points of the circular sector? "))
        print("")
        angleCircSecRad=np.radians(angleCircSecDeg)
        areaCircSec=0.5*(radiCircSec**2)*angleCircSecRad
        arclenCircSec=radiCircSec*angleCircSecRad
        print("The circular sector of radius ",radiCircSec," and angle (in degrees) of ",angleCircSecDeg," has an area of value ",areaCircSec," and an arc length of value ",arclenCircSec,".",sep="")
        circularsector=Circular_Sector(radius=radiCircSec,theta=angleCircSecDeg)
        circularsector.color="Black"
        circularsector.draw()
        break
    #Radius of a circle inscribed in a triangle
    if FigSelector==7:
        side1CircInsTri=abs(float(input("Length of side 1 of the triangle? ")))
        print("")
        side2CircInsTri=abs(float(input("Length of side 2 of the triangle? ")))
        print("")
        side3CircInsTri=abs(float(input("Length of side 3 of the triangle? ")))
        print("")
        angleCircInsTriDeg=float(input("Angle (in degrees) between side 1 and side 2 of the triangle? "))
        print("")
        angleCircInsTriRad=np.radians(angleCircInsTriDeg)
        semiperiCircInsTri=0.5*(side1CircInsTri+side2CircInsTri+side3CircInsTri)
        radiCircInsTri=(np.sqrt(semiperiCircInsTri*(semiperiCircInsTri-side1CircInsTri)*(semiperiCircInsTri-side2CircInsTri)*(semiperiCircInsTri-side3CircInsTri)))/(semiperiCircInsTri)
        x=sy.Symbol("x")
        y=sy.Symbol("y")
        eq1=x**2+y**2-side2CircInsTri**2
        eq2=(x-side1CircInsTri*np.cos(angleCircInsTriRad))**2+(y-side1CircInsTri*np.sin(angleCircInsTriRad))**2-side3CircInsTri**2
        point3=solve([eq1,eq2],[x,y],dict=True)
        values=list(point3[1].values())
        pointsCircInsTri=np.array([[0,0],[side1CircInsTri*np.cos(angleCircInsTriRad),side1CircInsTri*np.sin(angleCircInsTriRad)],[values[0],values[1]]])
        incenterxCircInsTri=(side1CircInsTri*values[0]+side2CircInsTri*side1CircInsTri*np.cos(angleCircInsTriRad)+side3CircInsTri*0)/(side1CircInsTri+side2CircInsTri+side3CircInsTri)
        incenteryCircInsTri=(side1CircInsTri*values[1]+side2CircInsTri*side1CircInsTri*np.sin(angleCircInsTriRad)+side3CircInsTri*0)/(side1CircInsTri+side2CircInsTri+side3CircInsTri)
        incenterCircInsTriCoord=np.array([[incenterxCircInsTri,incenteryCircInsTri]])
        print("The circle inscribed in a triangle of sides of length ",side1CircInsTri,", ",side2CircInsTri," and ",side3CircInsTri," and angle between the first two sides (in degrees) of ",angleCircInsTriDeg," has a semiperimeter of value ",semiperiCircInsTri," and the circle inscribed has a radius of value ",radiCircInsTri,".",sep="")
        circleinscribedintriangle=Circle_Inscribed_In_Triangle(points=pointsCircInsTri,radius=radiCircInsTri,incenterx=incenterxCircInsTri,incentery=incenteryCircInsTri)
        circleinscribedintriangle.color="Black"
        circleinscribedintriangle.draw()
        print("Vertex coordinates: ",pointsCircInsTri,sep="")
        print("Incenter coordinates: ",incenterCircInsTriCoord,sep="")
        break
    #Radius of a circle circumscribed by a triangle
    if FigSelector==8:
        side1CircCircuTri=abs(float(input("Length of side 1 of the triangle? ")))
        print("")
        side2CircCircuTri=abs(float(input("Length of side 2 of the triangle? ")))
        print("")
        side3CircCircuTri=abs(float(input("Length of side 3 of the triangle? ")))
        print("")
        angleCircCircuTriDeg=float(input("Angle (in degrees) between side 1 and side 2 of the triangle? "))
        print("")
        angleCircCircuTriRad=np.radians(angleCircCircuTriDeg)
        semiperiCircCircuTri=0.5*(side1CircCircuTri+side2CircCircuTri+side3CircCircuTri)
        radiCircCircuTri=(side1CircCircuTri*side2CircCircuTri*side3CircCircuTri)/(4*np.sqrt(semiperiCircCircuTri*(semiperiCircCircuTri-side1CircCircuTri)*(semiperiCircCircuTri-side2CircCircuTri)*(semiperiCircCircuTri-side3CircCircuTri)))
        x=sy.Symbol("x")
        y=sy.Symbol("y")
        eq1=x**2+y**2-side2CircCircuTri**2
        eq2=(x-side1CircCircuTri*np.cos(angleCircCircuTriRad))**2+(y-side1CircCircuTri*np.sin(angleCircCircuTriRad))**2-side3CircCircuTri**2
        point3=solve([eq1,eq2],[x,y],dict=True)
        values=list(point3[1].values())
        pointsCircCircuTri=np.array([[0,0],[side1CircCircuTri*np.cos(angleCircCircuTriRad),side1CircCircuTri*np.sin(angleCircCircuTriRad)],[values[0],values[1]]])
        xy1=pointsCircCircuTri[0]
        xy2=pointsCircCircuTri[1]
        xy3=pointsCircCircuTri[2]
        def circumcenter(points):
            (x1,y1),(x2,y2),(x3,y3)=points
            A=np.float64(np.array([[x3-x1,y3-y1],[x3-x2,y3-y2]]))
            Y=np.float64(np.array([(x3**2+y3**2-x1**2-y1**2),(x3**2+y3**2-x2**2-y2**2)]))
            if np.linalg.det(A)==0:
                return False
            Ainv=np.linalg.inv(A)
            X=0.5*np.dot(Ainv,Y)
            x,y=X[0],X[1]
            return (x,y)
        circumcenterCircCircuTriCoord=circumcenter((xy1,xy2,xy3))
        print("The circle circumscribed by a triangle of sides of length ",side1CircCircuTri,", ",side2CircCircuTri," and ",side3CircCircuTri," and angle between the first two sides (in degrees) of ",angleCircCircuTriDeg," has a semiperimeter of value ",semiperiCircCircuTri," and the circle circumscribed has a radius of value ",radiCircCircuTri,".",sep="")
        circlecircumscribedbytriangle=Circle_Circumscribed_By_Triangle(points=pointsCircCircuTri,radius=radiCircCircuTri,circumcenterx=circumcenterCircCircuTriCoord[0],circumcentery=circumcenterCircCircuTriCoord[1])
        circlecircumscribedbytriangle.color="Black"
        circlecircumscribedbytriangle.draw()
        print("Vertex coordinates: ",pointsCircCircuTri,sep="")
        print("Circumcenter coordinates: ",circumcenterCircCircuTriCoord,sep="")
        break
    #Area and perimeter of a regular polygon inscribed in a circle
    if FigSelector==9:
        NSidesPolInsCirc=abs(int(input("Number of sides of the regular polygon? ")))
        print("")
        radiPolInsCirc=abs(float(input("Radius of the circle? ")))
        print("")
        areaPolInsCirc=0.5*float(NSidesPolInsCirc)*(radiPolInsCirc**2)*np.sin((2*np.pi)/float(NSidesPolInsCirc))
        periPolInsCirc=2*float(NSidesPolInsCirc)*radiPolInsCirc*np.sin(np.pi/float(NSidesPolInsCirc))
        print("The regular polygon of ",NSidesPolInsCirc," sides, inscribed in a circle of radius ",radiPolInsCirc,", has an area of value ",areaPolInsCirc," and a perimeter of value ",periPolInsCirc,".",sep="")
        regpolygon=Regular_Polygon(radius=radiPolInsCirc,numVertices=NSidesPolInsCirc,circu=0)
        regpolygon.color="Black"
        regpolygon.draw()
        break
    #Area and perimeter of a regular polygon circumscribed by a circle
    if FigSelector==10:
        NSidesPolCircuCirc=abs(int(input("Number of sides of the regular polygon? ")))
        print("")
        radiPolCircuCirc=abs(float(input("Radius of the circle? ")))
        print("")
        areaPolCircuCirc=float(NSidesPolCircuCirc)*(radiPolCircuCirc**2)*np.tan(np.pi/float(NSidesPolCircuCirc))
        periPolCircuCirc=2*float(NSidesPolCircuCirc)*radiPolCircuCirc*np.tan(np.pi/float(NSidesPolCircuCirc))
        radiPolCircuCirc2=radiPolCircuCirc*float(mp.sec(np.pi/NSidesPolCircuCirc))
        print("The regular polygon of ",NSidesPolCircuCirc," sides, circumscribed by a circle of radius ",radiPolCircuCirc,", has an area of value ",areaPolCircuCirc," and a perimeter of value ",periPolCircuCirc,".",sep="")
        regpolygon=Regular_Polygon(radius=radiPolCircuCirc2,numVertices=NSidesPolCircuCirc,circu=1,radius2=radiPolCircuCirc)
        regpolygon.color="Black"
        regpolygon.draw()
        break
    #Area of a circle segment
    if FigSelector==11:
        radiCircSeg=abs(float(input("Radius of the circle? ")))
        print("")
        angleCircSegDeg=float(input("Central angle subtending the arc (in degrees) of the circle segment? "))
        print("")
        angleCircSegRad=np.radians(angleCircSegDeg)
        areaCircSeg=0.5*(radiCircSeg**2)*(angleCircSegRad-np.sin(angleCircSegRad))
        pointsTriCircSeg=np.array([[0,0],[radiCircSeg*np.cos((np.pi-angleCircSegRad)/2),radiCircSeg*np.sin((np.pi-angleCircSegRad)/2)],[-radiCircSeg*np.cos((np.pi-angleCircSegRad)/2),radiCircSeg*np.sin((np.pi-angleCircSegRad)/2)]])
        print("The circle segment of radius ",radiCircSeg," and central angle (in degrees) of ",angleCircSegDeg," has an area of value ",areaCircSeg,".",sep="")
        circlesegment=Circle_Segment(points=pointsTriCircSeg,radius=radiCircSeg,theta=angleCircSegDeg)
        circlesegment.color="Black"
        circlesegment.draw()
        break
    #Area and arc length of a parabolic segment
    if FigSelector==12:
        maxheightParSeg=abs(float(input("Maximum height of the parabolic segment? ")))
        print("")
        baseParSeg=abs(float(input("Length of the base of the parabolic segment? ")))
        print("")
        areaParSeg=(2/3)*maxheightParSeg*baseParSeg
        arclenParSeg=0.5*np.sqrt(baseParSeg**2+16*(maxheightParSeg)**2)+((baseParSeg**2)/(8*maxheightParSeg))*np.log((4*maxheightParSeg+np.sqrt(baseParSeg**2+16*(maxheightParSeg)**2))/(baseParSeg))
        x_lista_ParSeg=np.linspace(-baseParSeg/2,baseParSeg/2,100)
        y_lista_ParSeg=[]
        for x in x_lista_ParSeg:
            y=maxheightParSeg*(1-((x**2)/((baseParSeg/2)**2)))
            y_lista_ParSeg.append(y)
        print("The parabolic segment of maximum height of ",maxheightParSeg," and base of length ",baseParSeg," has an area of value ",areaParSeg," and an arc length of value ",arclenParSeg,".",sep="")
        plt.plot(x_lista_ParSeg,y_lista_ParSeg,marker="o",markevery=[0,y_lista_ParSeg.index(max(y_lista_ParSeg)),99],color="Black")
        plt.title("Parabolic_Segment")
        plt.grid()
        plt.tight_layout()
        plt.show()
        break
    #Area and perimeter of an ellipse
    if FigSelector==13:
        semihoraxisEll=abs(float(input("Length of the semi-horitzontal axis of the ellipse? ")))
        print("")
        semiveraxisEll=abs(float(input("Length of the semi-vertical axis of the ellipse? ")))
        print("")
        if semihoraxisEll==semiveraxisEll:
            print("It seems you want to calculate the area and perimeter of a *circle*, let me help you with that...")
            print("")
            radiCirc=semihoraxisEll
            areaCirc=np.pi*radiCirc**2
            periCirc=2*np.pi*radiCirc
            print("The *circle* of radius ",radiCirc," has an area of value ",areaCirc," and a perimeter of value ",periCirc,".",sep="")
            circle=Circle(radius=radiCirc)
            circle.color="Black"
            circle.draw()
            break
        areaEll=np.pi*semihoraxisEll*semiveraxisEll
        theta=sy.Symbol("theta")
        periEll=4*semihoraxisEll*sy.integrate(sy.sqrt(1-(((sy.sqrt(semihoraxisEll**2-semiveraxisEll**2))/(semihoraxisEll))**2)*(sy.sin(theta))**2),(theta,0,np.pi/2))
        print("The ellipse of semi-horitzontal axis of length ",semihoraxisEll," and semi-vertical axis of length ",semiveraxisEll," has an area of value ",areaEll," and a perimeter of value ",periEll,".",sep="")
        ellipse=Ellipse(width=2*semihoraxisEll,height=2*semiveraxisEll)
        ellipse.color="Black"
        ellipse.draw()
        break