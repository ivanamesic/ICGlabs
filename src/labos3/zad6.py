from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse
import math
import numpy
import random
import sys

class Point:
    def __init__(self, x, y, z):
        self.x=x
        self.y=y
        self.z=z
class Coefs:
    def __init__(self,a,b,c,d):
        self.a=a
        self.b=b
        self.c=c
        self.d=d
class Vertex:
    def __init__(self):
        self.sumx=self.sumy=self.sumz=0
        self.n=0
        self.Nx=self.Ny=self.Nz=0
window = pyglet.window.Window()
width=window.get_size()[0]
height=window.get_size()[1]
G_default=Point(0.0,0.,0.)
O_default=Point(0,-9,26.)
O=O_default
G=G_default

vertices=[]
polygons=[]
trajectory=[]
coefs=[]
v_details=[]
polygon_centers=[]
t=0
step=0.1

def readInput(file):
    with open(file, 'r') as fp:
        for line in fp:
            if (line.startswith("#") or len(line)==0):continue
            if (line.startswith("v")):
                line=line[2:]
                vertices.append([float(j) for j in line.split(" ")])
            if (line.startswith("f")):
                line=line[2:]
                polygons.append([int(j)-1 for j in line.split(" ")])
    convert_coordinates(vertices)


def read_trajectory(file):
    with open(file, 'r') as fp:
        for line in fp:
            trajectory.append([float(j) for j in line.split(" ")])
    convert_coordinates(trajectory)

@window.event
def on_show():
    glClearColor(1.0, 1.0, 1.0, 1.0)
    window.set_caption("Zadatak 6: Bezierove krivulje")


@window.event
def on_draw():
    glViewport(0,0,width,height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    gluOrtho2D(-1.5,1.5,-1.5,1.5)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

    glClearColor(0.,0.,0., 1.0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glEnable(GL_DEPTH_TEST)
    glPointSize(1.0)

    draw_polygons()
    draw_bezier_line(trajectory)

    #draw_with_animation()

@window.event
def on_resize(w,h):
    width=w
    height=h
    glViewport(0,0,width,height)
    update_scene()

def update_scene():
    global O
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, width/height, 0.5, 8.0)
    gluLookAt (O.x, O.y, O.z, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    glMatrixMode (GL_MODELVIEW)

def draw_polygons():
    for i in range (len(polygons)):
        c=coefs[i]
        p=polygons[i]

        cos = c.a * O.x + c.b * O.y + c.c * O.z + c.d

        if (cos > 0):
            #glColor3f(random.random(),random.random(),random.random())
            glColor3f(1.0,1.0,1.0 )
            glBegin(GL_LINE_LOOP)
            for i in p:
                v = vertices[i]
                T = numpy.array([v[0], v[1], v[2], 1])
                A_s = T.dot(T1_5).dot(P)
                glVertex2f(A_s.item(0, 0), A_s.item(0, 1))
            glEnd()

def center(pol):
    v1 = vertices[pol[0]]
    v2 = vertices[pol[1]]
    v3 = vertices[pol[2]]
    return [(v1[0]+v2[0]+v3[0])/3, (v1[1]+v2[1]+v3[1])/3, (v1[2]+v2[2]+v3[2])/3]

t = 0;step=0.1


def update_glediste(dt):
    global t,step,G, O
    t += step
    if t >= 1: step = -0.1
    if t <= 0: step = 0.1
    n = len(trajectory) - 1
    p = numpy.array([0, 0, 0])
    for i in range(n + 1):
        b = bernstein_function(n, i, t)
        res = numpy.array(trajectory[i]) * b
        p = numpy.add(numpy.array(p), res)
    G= Point(p[0], p[1], p[2])
    view_transformation()

#Racunanje i iscrtavanje krivulje Beziera definirane Bersteinovim tezinskim polinomima
def draw_bezier_line(polygon):
    glColor3f(0.0, 0.0, 1.0 )
    glBegin(GL_LINE_LOOP)
    for i in polygon:
        glVertex3f(i[0],i[1],i[2])
    glEnd()
    n = len(polygon)-1
    global O
    glColor3f(1.0,0.0,0.0)
    glBegin(GL_LINE_STRIP)
    for t in numpy.arange(0.,1.,0.01):
        p=numpy.array([0,0,0])
        for i in range(n+1):
            b=bernstein_function(n,i,t)
            res=numpy.array(polygon[i])*b
            p=numpy.add(numpy.array(p),res)
        glVertex3f(p[0],p[1],p[2])
    glEnd()


#racunanje bazne funkcije Bernsteinovim polinomima
def bernstein_function(n,i,t):
    x=math.factorial(n)/(math.factorial(i)*math.factorial(n-i))
    return x*math.pow(t,i)*math.pow(1-t, n-i)

#pretvaranje koordinata (pomak i skaliranje)
def convert_coordinates(array):
    max_x=min_x=array[0][0]
    max_y=min_y=array[0][1]
    max_z=min_z=array[0][2]
    for i in array:
        max_x=max(max_x,i[0])
        min_x=min(min_x,i[0])
        max_y=max(max_y,i[1])
        min_y=min(min_y,i[1])
        max_z=max(max_z,i[2])
        min_z=min(min_z,i[2])

    xs = (max_x + min_x) / 2.0
    ys = (max_y + min_y) / 2.0
    zs = (max_z + min_z) / 2.0

    m = max(max_x - min_x, max_y - min_y, max_z - min_z)
    global O, G
    k = 2 / m

    for v in array:
        v[0]=(v[0]-xs)*k
        v[1]=(v[1]-ys)*k
        v[2]=(v[2]-zs)*k
    O.x = (O.x - xs) * k
    O.y = (O.y - ys) * k
    O.z = (O.z - zs) * k
    G.x = (G.x - xs) * k
    G.y = (G.y - ys) * k
    G.z = (G.z - zs) * k


#racunanje koeficijenata bridova
def calculate_coefs():
    for pol in polygons:
        v1 = vertices[pol[0]]
        v2 = vertices[pol[1]]
        v3 = vertices[pol[2]]
        a = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1])
        b = -(v2[0] - v1[0]) * (v3[2] - v1[2]) + (v2[2] - v1[2]) * (v3[0] - v1[0])
        c = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])
        d = -v1[0] * a - v1[1] * b - v1[2] * c
        coefs.append(Coefs(a,b,c,d))
        polygon_centers.append(center(pol))

def view_transformation():
    T1=numpy.matrix([[1,0,0,0],[0,1,0,0], [0,0,1,0],[-O.x,-O.y,-O.z,1]])
    G1=Point(G.x-O.x,G.y-O.y,G.z-O.z)

    sin_alfa=G1.y/math.sqrt(G1.x**2+G1.y**2)
    cos_alfa=G1.x/math.sqrt(G1.x**2+G1.y**2)
    T2=numpy.matrix([[cos_alfa,-sin_alfa,0,0],[sin_alfa,cos_alfa,0,0], [0,0,1,0],[0,0,0,1]])

    G2=Point(math.sqrt(G1.x**2+G1.y**2),0,G1.z)

    sin_beta=G2.x/math.sqrt(G2.x**2+G2.z**2)
    cos_beta=G2.z/math.sqrt(G2.x**2+G2.z**2)
    T3=numpy.matrix([[cos_beta,0,sin_beta,0],[0,1,0,0], [-sin_beta,0,cos_beta,0],[0,0,0,1]])

    G3=Point(0,0,math.sqrt(G2.x**2+G2.z**2))

    T4=numpy.matrix([[0,-1,0,0], [1,0,0,0],[0,0,1,0], [0,0,0,1]])
    T5=numpy.matrix([[-1,0,0,0], [0,1,0,0],[0,0,1,0], [0,0,0,1]])
    T2_5=T2.dot(T3).dot(T4).dot(T5)
    global T1_5
    T1_5=T1.dot(T2_5)
    perspective_projection()

def perspective_projection():
    H=math.sqrt(math.pow(O.x-G.x,2)+math.pow(O.y-G.y,2)+math.pow(O.z-G.z,2))
    global P,k
    k=1/H
    P=numpy.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,k,0]])

#citanje ulaznih datoteka i racunanje pocetnih parametara
putanja_file=sys.argv[2]
vrhovi_file=sys.argv[1]
read_trajectory(putanja_file)
readInput(vrhovi_file)
calculate_coefs()
view_transformation()
pyglet.clock.schedule_interval(update_glediste, 1 / 10.0)
pyglet.app.run();