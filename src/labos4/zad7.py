from pyglet.gl import *
import sys
import numpy
import math

class Coefs:
    def __init__(self,a,b,c,d):
        self.a=a
        self.b=b
        self.c=c
        self.d=d
class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
class Polygon:
    def __init__(self, v):
        self.vertices = v
    def set_intensity(self, i):
        self.intensity = i
    def set_coefs(self, coef, norm):
        self.coef = coef
        self.norm = norm
    def set_isHidden(self, isHidden):
        self.isHidden=isHidden
    def set_center(self, center):
        self.center = center

window = pyglet.window.Window()
width = window.get_size()[0]
height = window.get_size()[1]
G=G_default = Point(0,0,0)
O=O_default = Point(6,2,6)
c_shading = True
wired = False

vertices = []
polygons = []
vertexToPolygonMap = {}
vertexNormals = {}

# Postavljanje očišta i gledišta na default vrijednosti
def define_scene():
    global O, G
    O = Point(O_default.x, O_default.y, O_default.z)
    G = Point(G_default.x, G_default.y, G_default.z)



@window.event
def on_resize(w, h):
    width = w
    height = h
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluOrtho2D(-1, 1, -1, 1)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glViewport(0, 0, width, height)


@window.event
def on_draw():
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluOrtho2D(-2, 2, -2, 2)

    glClearColor(1.0, 1.0, 1.0, 1.0)
    glClear(GL_COLOR_BUFFER_BIT)
    glPointSize(1.0)

    if wired:
        draw_body_wired()
    else:
        draw_body_with_colors()

#--------------------------------------------------------------------------------------
# 0. Učitavanje podataka iz datoteke
def read_file(file):
    with open(file, 'r') as fp:
        for line in fp:
            if (line.startswith("#") or len(line) == 0): continue
            if (line.startswith("v")):
                line = line[2:]
                vertices.append([float(j) for j in line.split(" ")])
            if (line.startswith("f")):
                line = line[2:]
                v = [int(j) - 1 for j in line.split(" ")]
                polygons.append(Polygon(v))

#--------------------------------------------------------------------------------------
# 1a. Konverzija koordinata vrhova (pomak i skaliranje)

def convert_coordinates():
    global max_x, max_y, max_z, min_x, min_y, min_z
    max_x = min_x = vertices[0][0]
    max_y = min_y = vertices[0][1]
    max_z = min_z = vertices[0][2]
    for i in vertices:
        max_x = max(max_x, i[0])
        min_x = min(min_x, i[0])
        max_y = max(max_y, i[1])
        min_y = min(min_y, i[1])
        max_z = max(max_z, i[2])
        min_z = min(min_z, i[2])

    xs = (max_x + min_x) / 2.0
    ys = (max_y + min_y) / 2.0
    zs = (max_z + min_z) / 2.0

    m = max(max_x - min_x, max_y - min_y, max_z - min_z)
    k = 2 / m
    max_x = (max_x - xs) * k
    max_y = (max_y - ys) * k
    max_z = (max_z - zs) * k
    min_x = (min_x - xs) * k
    min_y = (min_y - ys) * k
    min_z = (min_z - zs) * k

    for v in vertices:
        v[0] = (v[0] - xs) * k
        v[1] = (v[1] - ys) * k
        v[2] = (v[2] - zs) * k

#--------------------------------------------------------------------------------------
# 1b. Racunanje koeficijenata bridova i norme poligona
def calculate_coefs():
    for p in polygons:
        pol = p.vertices
        v1 = vertices[pol[0]]
        v2 = vertices[pol[1]]
        v3 = vertices[pol[2]]
        xs = (v1[0]+v2[0]+v3[0])/3
        ys = (v1[1]+v2[1]+v3[1])/3
        zs = (v1[2]+v2[2]+v3[2])/3
        p.set_center(Point(xs, ys, zs))
        a = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1])
        b = -(v2[0] - v1[0]) * (v3[2] - v1[2]) + (v2[2] - v1[2]) * (v3[0] - v1[0])
        c = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])
        d = -v1[0] * a - v1[1] * b - v1[2] * c
        norm = math.sqrt(a**2+b**2+c**2)
        normVector=[a/norm, b/norm, c/norm]
        p.set_coefs(Coefs(a,b,c,d),normVector)

#--------------------------------------------------------------------------------------
# 1c. Transformacija pogleda i perspektivna projekcija
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
    global P
    k=1/H
    P=numpy.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,k,0]])

@window.event
def on_key_press(symbol, modifiers):
    global c_shading, Iz, O, wired
    if symbol == pyglet.window.key.R:
        print("Postavi ociste na nulu")
        O.x = 0.0
    elif symbol == pyglet.window.key.O:
        print("Povećavam x koordinatu očišta")
        O.x = O.x + 0.1
    elif symbol == pyglet.window.key.P:
        print("Smanjujem x koordinatu očišta")
        O.x = O.x - 0.1
    elif symbol == pyglet.window.key.L:
        print("Povećavam x koordinatu izvora")
        Iz.x = Iz.x + 0.1
    elif symbol == pyglet.window.key.K:
        print("Smanjujem x koordinatu izvora")
        Iz.x = Iz.x - 0.1
    elif symbol == pyglet.window.key.U:
        Iz.y = Iz.y + 0.1
    elif symbol == pyglet.window.key.I:
        Iz.y = Iz.y - 0.1
    elif symbol == pyglet.window.key.C:
        print("Mijenjam tip sjencanja")
        if c_shading:
            c_shading=False
        else:
            c_shading=True
    elif symbol == pyglet.window.key.Z:
        print("Mijenjam tip crtanja tijela")
        if wired:
            wired=False
        else:
            wired=True
    remove_polygons()
    if c_shading:
        calculate_constant_intensities()
    else:
        calculate_gouraud_intensity()
    view_transformation()
#--------------------------------------------------------------------------------------
# 2. Odredivanje straznjih poligona
def remove_polygons():
    for pol in polygons:
        p = pol.vertices
        c = pol.coef
        cos = c.a * O.x + c.b * O.y + c.c * O.z + c.d
        if (cos > 0):
            pol.set_isHidden(False)
        else:
            pol.set_isHidden(True)
        for j in range(3):
            if (p[j] in vertexToPolygonMap.keys()):
                vertexToPolygonMap.get(p[j]).append(pol)
            else:
                vertexToPolygonMap[p[j]] = [pol]

#--------------------------------------------------------------------------------------
# 3. Iscrtavanje tijela uz micanje straznjih poligona
def draw_body_with_colors():
    for p in polygons:
        if (p.isHidden): continue
        glBegin(GL_TRIANGLES)
        for i in range(len(p.vertices)):
            ind = p.vertices[i]
            v = vertices[ind]
            T = numpy.array([v[0], v[1], v[2], 1])
            A_s = T.dot(T1_5).dot(P)
            I = (int)(p.intensity[i])
            glColor3ub(I,0,0)
            glVertex2f(A_s.item(0, 0), A_s.item(0, 1))
        glEnd()

def draw_body_wired():
    for p in polygons:
        if (p.isHidden): continue
        glColor3f(0.0, 0.0, 0.0)
        glBegin(GL_LINE_LOOP)
        for i in range(len(p.vertices)):
            ind = p.vertices[i]
            v = vertices[ind]
            T = numpy.array([v[0], v[1], v[2], 1])
            A_s = T.dot(T1_5).dot(P)
            glVertex2f(A_s.item(0, 0), A_s.item(0, 1))
        glEnd()
#--------------------------------------------------------------------------------------
# 4. Zadavanje polozaja Izvora
Iz = Point(6., 6., 10.)

#--------------------------------------------------------------------------------------
# 5. Odredivanje intenziteta poligona

Ia = 80
ka = 0.6
Ii = 220
kd = 0.7

def calculate_constant_intensities():
    for p in polygons:
        center=p.center
        L=numpy.array([Iz.x-center.x, Iz.y-center.y, Iz.z-center.z])
        LN = numpy.inner(L,p.norm)
        if LN < 0:
            dif = 0
        else:
            dif = Ii*kd*LN
        I = Ia * ka + dif;
        p.set_intensity([I,I,I])

#--------------------------------------------------------------------------------------
# 7. Racunanje normala vrhova
def calculate_vertex_normals():
    for v in vertexToPolygonMap.keys():
        n = len(vertexToPolygonMap.get(v))
        res = numpy.array([0.0,0.0,0.0])
        for n_i in vertexToPolygonMap.get(v):
            res = res + numpy.array(n_i.norm)
        norm = numpy.linalg.norm(res)
        res = (res/norm )/ n
        vertexNormals[v]=res

#--------------------------------------------------------------------------------------
# 8. Racunanje intenziteta vrhova
def calculate_gouraud_intensity():
    for p in polygons:
        inten = []
        for v in p.vertices:
            vertex = vertices[v]
            L = numpy.array([Iz.x-vertex[0], Iz.y-vertex[1], Iz.z-vertex[2]])
            LN = numpy.inner(L,vertexNormals.get(v))
            if LN < 0:
                dif = 0
            else:
                dif = Ii*kd*LN
            inten.append(Ia * ka + dif)
        p.set_intensity(inten)


# Ucitavanje datoteke, izračun početnih postavki
file = sys.argv[1]
read_file(file)
convert_coordinates()
calculate_coefs()
view_transformation()
remove_polygons()
calculate_vertex_normals()
calculate_constant_intensities()
pyglet.app.run()