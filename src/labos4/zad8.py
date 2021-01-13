from pyglet.gl import *
import cmath

print("eps=100;m=16;u_min=-1.5;u_max=0.5;v_min=-1;v_max=1")

eps=int(input("Upišite prag(eps): "))
m=int(input("Upišite maksimalni broj iteracija: "))
u_min=float(input("u_min: "))
u_max=float(input("u_max: "))
v_min=float(input("v_min: "))
v_max=float(input("v_max: "))
window=pyglet.window.Window()
x_max,y_max=[400,400]

@window.event
def on_show():
    window.set_size(x_max,y_max)
    glClearColor(1.0, 1.0, 1.0, 1.0)
    window.set_caption("Zadatak 8.3.1: Mandelbrotov fraktal")

@window.event
def on_resize(w,h):
    x_max=w
    y_max=h
    glMatrixMode(gl.GL_PROJECTION)
    gluOrtho2D(0,x_max-1,0,y_max-1)
    glLoadIdentity()
    glMatrixMode(gl.GL_MODELVIEW)
    glViewport(0,0,x_max,y_max)

@window.event
def on_draw():
    glClearColor(1.0, 1.0, 1.0, 1.0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glPointSize(1.0)
    glColor3f(1.,0.,0.)
    glBegin(GL_LINE_LOOP)
    glVertex2i(0,100)
    glVertex2i(100,0)
    glEnd()
    draw_mandelbrot()

def draw_mandelbrot():
    glBegin(GL_POINTS)
    for x in range(x_max):
        for y in range(y_max):
            u0=(u_max-u_min)*x/x_max + u_min
            v0=(v_max-v_min)*y/y_max + v_min
            k=-1;c_real=u0;c_imag=v0;z=complex(0,0)
            c=complex(c_real,c_imag)
            r = abs(z)
            while(r < eps and k < m):
                k += 1
                z=z**2+c
                r=abs(z)
            colorPoint(k)
            glVertex2i(x,y)
    glEnd()

def colorPoint(k):
    if k == m:
        glColor3f(0., 0., 0.)
    else:
        if m < 32:
            limit = m
        else:
            limit = 32
        glColor3f(k / limit, 1 - k / limit / 2., 0.8 - k / limit / 3.)

pyglet.app.run()
