from pyglet.gl import *
from pyglet.window import key
from pyglet.window import mouse
import numpy
import math
import sys

class Point:
    def __init__(self, x, y, z):
        self.x=x
        self.y=y
        self.z=z
        
window = pyglet.window.Window()
width=window.get_size()[0]
height=window.get_size()[1]
G_default=Point(0.0,0.,0.)
O_default=Point(1.,1, 1)

vertices=[]
polygons=[]

global G,O
@window.event
def on_show():
    print("Koordinate gledišta: (%f, %f, %f); koordinate očišta: (%f, %f, %f)" % (G.x, G.x, G.z, O.x, O.y, O.z))
    print("Za povećanje/smanjenje x koordinate očišta: pritisnuti O/shift+O")
    print("Za povećanje/smanjenje y koordinate očišta: pritisnuti P/shift+P")
    print("Za povećanje/smanjenje z koordinate očišta: pritisnuti Q/shift+Q")
    print("Za povećanje/smanjenje x koordinate gledišta: pritisnuti X/shift+X")
    print("Za povećanje/smanjenje y koordinate gledišta: pritisnuti Y/shift+Y")
    print("Za povećanje/smanjenje z koordinate gledišta: pritisnuti Z/shift+Z")
    print("Za vraćanje na početne postavke stisnuti Backspace")
    print("Za ispis koordinata očišta i gledišta, stisnuti S")

#Postavljanje očišta i gledišta na default vrijednosti
def define_scene():
    global O,G
    O=Point(O_default.x,O_default.y,O_default.z)
    G=Point(G_default.x,G_default.y,G_default.z)

@window.event
def on_key_press(symbol, modifiers):
    if symbol == key.BACKSPACE:
        define_scene()
        print("Vraćam gledište i očište na početne vrijednosti...")
    elif symbol == key.S:
        print("Očište: (%f,%f,%f)" % (O.x,O.y,O.z))
        print("Gledište: (%f,%f,%f)" % (G.x,G.y,G.z))
    elif symbol == key.O:
        if modifiers == key.MOD_SHIFT:
            s="Smanjujem x koordinatu očišta"
            r = O.x - 0.1
        else:
            s="Povećavam x koordinatu očišta"
            r = O.x + 0.1
        if (is_in_shape(r,O.y,O.z)):
            print("Očište ne smije biti unutar objekta, neće se dogoditi promjena")
        else:
            print(s)
            O.x = r
    elif symbol == key.P:
        if modifiers == key.MOD_SHIFT:
            s="Smanjujem y koordinatu očišta"
            r = O.y - 0.1
        else:
            s="Povećavam y koordinatu očišta"
            r = O.y + 0.1
        if (is_in_shape(O.x, r, O.z)):
            s="Očište ne smije biti unutar objekta, neće se dogoditi promjena"
        else:
            O.y = r
            print(s)
    elif symbol == key.Q:
        if modifiers == key.MOD_SHIFT:
            s="Smanjujem z koordinatu očišta"
            r = O.z - 0.5
        else:
            s="Povećavam z koordinatu očišta"
            r = O.z + 0.5
        if (is_in_shape(O.x, O.y, r)):
            print("Očište ne smije biti unutar objekta, neće se dogoditi promjena")
        else:
            O.z = r
            print(s)
    elif symbol == key.X:
        if modifiers == key.MOD_SHIFT:
            print("Smanjujem x koordinatu gledišta...")
            G.x -= 0.1
        else:
            print("Povećavam x koordinatu gledišta...")
            G.x += 0.1
    elif symbol == key.Y:
        if modifiers == key.MOD_SHIFT:
            print("Smanjujem y koordinatu gledišta...")
            G.y -=0.1
        else:
            print("Povećavam y koordinatu gledišta...")
            G.y += 0.1
    elif symbol == key.Z:
        if modifiers == key.MOD_SHIFT:
            print("Smanjujem z koordinatu gledišta...")
            G.z -=0.1
        else:
            print("Povećavam z koordinatu gledišta...")
            G.z += 0.1
    view_transformation()

#Provjera nalazi li se tocka unutar lika
def is_in_shape(x,y,z):
    global max_x,max_y,max_z,min_x,min_y,min_z
    return (x>min_x and x < max_x and y > min_y and y < max_y and z > min_z and z < max_z)

@window.event
def on_draw():
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    gluOrtho2D(-2, 2, -2, 2)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glViewport(0, 0, width, height)

    glClearColor(1.0, 1.0, 1.0, 1.0)
    glClear(GL_COLOR_BUFFER_BIT)
    glPointSize(1.0)

    draw_polygons_rotated()

#Iscrtavanje poligona
def draw_polygons_rotated():
    glColor3f(0.0,0.0,0.0)
    for p in polygons:
        glBegin(GL_LINE_LOOP)
        for i in p:
            v=vertices[i]
            T=numpy.array([v[0],v[1],v[2],1])
            A_s=T.dot(T1_5).dot(P)
            glVertex2f(A_s.item(0,0),A_s.item(0,1))
        glEnd()

#Transformacija pogleda
def view_transformation():
    global G
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

#Računanje matrice perspektivne projekcije
def perspective_projection():
    H=math.sqrt(math.pow(O.x-G.x,2)+math.pow(O.y-G.y,2)+math.pow(O.z-G.z,2))
    global P
    k=1/H
    P=numpy.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,k,0]])

@window.event
def on_resize(w,h):
    width=w
    height=h
    glViewport(0,0,width,height)
    glMatrixMode(gl.GL_PROJECTION)
    gluOrtho2D(-2,2,-2,2)
    glLoadIdentity()
    glMatrixMode(gl.GL_MODELVIEW)

#Učitavanja podataka iz datoteke
def read_file(file):
    with open(file, 'r') as fp:
        for line in fp:
            if (line.startswith("#") or len(line)==0):continue
            if (line.startswith("v")):
                line=line[2:]
                vertices.append([float(j) for j in line.split(" ")])
            if (line.startswith("f")):
                line=line[2:]
                polygons.append([int(j)-1 for j in line.split(" ")])

#Konverzija koordinata vrhova (pomak i skaliranje)
def convert_coordinates():
    global max_x,max_y,max_z,min_x,min_y,min_z
    max_x=min_x=vertices[0][0]
    max_y=min_y=vertices[0][1]
    max_z=min_z=vertices[0][2]
    for i in vertices:
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
    k = 2 / m
    max_x=(max_x-xs)*k
    max_y=(max_y-ys)*k
    max_z=(max_z-zs)*k
    min_x=(min_x-xs)*k
    min_y=(min_y-ys)*k
    min_z=(min_z-zs)*k

    for v in vertices:
        v[0]=(v[0]-xs)*k
        v[1]=(v[1]-ys)*k
        v[2]=(v[2]-zs)*k
    global O, G
    O.x = (O.x - xs) * k
    O.y = (O.y - ys) * k
    O.z = (O.z - zs) * k
    G.x = (G.x - xs) * k
    G.y = (G.y - ys) * k
    G.z = (G.z - zs) * k
#Ucitavanje datoteke, izračun početnih postavki
define_scene()
file = sys.argv[1]
read_file(file)
convert_coordinates()
view_transformation()
pyglet.app.run()
