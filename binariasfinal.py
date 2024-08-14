import numpy as np
import kepler
from vpython import *
import os


#Creamos una función para leer el dato ingresado por el usuario
# y proporcionar un valor predeterminado si la entrada está vacía.
def obtener_entrada(instrucciones, valor_predeterminado):
    entrada = input(instrucciones)
    return float(entrada) if entrada.strip() else valor_predeterminado

#Ingresamos los datos:
M1 = obtener_entrada('Ingrese la masa en masas solares (valor predeterminado: 3.0)', 3.0) 
while M1<0: M1 = obtener_entrada('La masa debe ser mayor a cero, ingrese la masa en masas solares nuevamente ', 3.0)
M2 = obtener_entrada('Ingrese la segunda masa en masas solares (valor predeterminado: 1.0)', 1.0)
while M2<0: M2 = obtener_entrada('La masa debe ser mayor a cero, ingrese la masa en masas solares nuevamente ', 1.0)
i = obtener_entrada('Ingrese la inclinación en grados (valor predeterminado: 90.0)', 90.0) * np.pi / 180  # Convertir a radianes
P = obtener_entrada('Ingrese el período en años (valor predeterminado: 0.019)', 0.019)
while P<0: P = obtener_entrada('El periodo debe ser positivo, ingrese el período en años nuevamente ', 0.019)
e = obtener_entrada('Ingrese la excentricidad (valor predeterminado: 0.5)', 0.5)
while e<=0 or e>=1: e = obtener_entrada('La excentricidad debe estar entre 0 y 1, ingrese la excentricidad nuevamente ', 0.5)
w1 = obtener_entrada('Ingrese el argumento del periastro en grados (valor predeterminado: 229.0)', 229.0) * np.pi / 180  # Convertir a radianes
V0 = obtener_entrada('Ingrese la velocidad baricentral en km/s (valor predeterminado: -6.91)', -6.91) *1000 #Convertir a m/s



#Calculemos el semieje mayor de la orbita relativa:
a=((P**2*(M1+M2))**(1/3))*150e9

#Calculemos la longitud del periastro para la secundaria
w2=w1-np.pi

#Calculamos los semiejes mayores de cada orbita
a2=a*M1/(M1+M2)
a1=a-a2

#Calculemos la posición inicial de la primaria
r1=a1*(1-e)

#Calculemos la posición inicial de la secundaria
r2=(r1/(M2/M1))

#Pasamos el periodo a segundos para los siguientes cálculos
Pseg=P*3.154e7


#Calculamos las semiamplitudes de cada componente
k1=(2*np.pi*a1*np.sin(i))/(Pseg*np.sqrt(1-e**2))
k2=(2*np.pi*a2*np.sin(i))/(Pseg*np.sqrt(1-e**2))


#Definimos la ecuación de las tangentes 
def eq_tang(E):
    return 2*np.arctan2(np.sqrt((1+e)/(1-e))*np.tan(E/2.),1)
 

#Definimos la ecuación de velocidad radial de cada componente
def eq_Vr1(var):
    return V0+k1*(e*np.cos(w1)+np.cos(w1+var))

def eq_Vr2(var):
    return V0+k2*(e*np.cos(w2)+np.cos(w2+var))


#Animación orbitas

#Creamos un lienzo en el cuál irán las órbitas
c1=canvas(title='Sistema binario',width=700,height=500,align='left')
G = 6.7e-11 # constante de gravitación universal

#Calculamos el momento angular específico de la primaria
L1=sqrt(4*pi**2*a1**4*(1-e**2)/(Pseg**2))

#Reemplazamos en la siguiente ecuación para obtener la velocidad inicial de la primaria
v1=sqrt((L1**2/(a1*(1-e**2)))*(2/r1-1/a1))
print('la velocidad inicial de 1 es',v1)

#Hacemos lo mismo para la secundaria
L2=sqrt(4*pi**2*a2**4*(1-e**2)/(Pseg**2))
v2=sqrt((L2**2/(a2*(1-e**2)))*(2/r2-1/a2))
print('la velocidad inicial de 2 es',v2)

#Creamos una función para definir el radio según la masa en metros
def radioSegunMasa(m:float):
    if m<1.66:
        rad=1.06*m**(0.945)
    else :
        rad=1.33*m**(0.555)   
    return rad*6.957e8



# Función para ajustar el retain según el tamaño de la órbita
def calcular_retain(tamano_orbita):
    base_retain = 1000  # Base del retain para órbitas pequeñas
    factor_retain = 5000  # Factor de escala para órbitas grandes
    retain_value = int(base_retain + factor_retain * tamano_orbita)
    return retain_value

# Ajustar el valor de retain según el tamaño de la órbita
tamano_orbita = a / 1e11  # Escalar el tamaño de la órbita para ajustar retain
retain_value = calcular_retain(tamano_orbita)



#Definimos la estrella primaria
est1 = sphere(canvas=c1,pos=vector(r1,0,0), radius=radioSegunMasa(M1), color=color.magenta, 
                make_trail=True,  interval=5, retain=retain_value)
est1.masa = M1*1.988e30 #Pasamos la masa a kg
est1.p = vector(0, v1, 0) * est1.masa #Momento lineal primaria

#Definimos la estrella secundaria
est2 = sphere(canvas=c1,pos=vector(-r2,0,0), radius=radioSegunMasa(M2), color=color.green,
                make_trail=True, interval=5, retain=retain_value)
est2.masa = M2*1.988e30 #Pasamos la masa a kg
est2.p = vector(0, -v2, 0) * est2.masa #Momento lineal secundaria

#Agregamos las correspondientes descripciones 
label1=label(pos=vec(r1, 0, 0), text='primaria'if M1>M2 else 'secundaria', xoffset=20, yoffset=50, space=30, 
    height=16, border=4, font='sans')
label2=label(pos=vec(-r2, 0, 0), text='secundaria' if M1>M2 else 'primaria', xoffset=20, yoffset=50, space=30,
    height=16, border=4, font='sans')


#Calculamos las energías:

#Energía cinética:
Ec =  0.5*est1.masa*v1**2+0.5*est2.masa*v2**2
print('La energía cinetica es',Ec)

#Energía potencial:
Ep = - G*est1.masa*est2.masa/(mag(est1.pos-est2.pos))
print('La energía potencial es',Ep)

#Energía total:
E=Ec+Ep
print('La energía total es',E)

#Tenemos en cuenta que
if E>0:
    print('El sistema no está ligado, elija otros parámetros')
    quit()

#Momento angular:
L= cross(est1.pos,est1.p)+cross(est2.pos,est2.p)
print('el momento angular es',L)



#Definimos las curvas de velocidad radial para cada componente

vr1_grafico = graph(align='left',title='Curvas de velocidad radial',
                   xtitle="Fase", ytitle = "Vr [m/s]", width=700,
                     height = 250, xmin = 0, xmax = 10)
vr2_grafico = graph(align='left',
                   xtitle="Fase", ytitle = "Vr [m/s]", width=700,
                     height = 250, xmin = 0, xmax = 10)
vr1_curva = gcurve(color=color.magenta,graph=vr1_grafico,
                   label='Estrella primaria' if M1>M2 else 'Estrella secundaria')
vr2_curva = gcurve(color=color.green,graph=vr2_grafico,
                   label='Estrella secundaria' if M1>M2 else 'Estrella primaria')


#Creamos un botón para pausar la animación
corriendo = True
def Pausa(b): # b = botón
    global corriendo, recordar_dt, dt
    corriendo = not corriendo
    if corriendo:
        b.text = "Pausa"
        dt = recordar_dt
    else: 
        b.text = "Iniciar"
        recordar_dt = dt
        dt = 0
    return

button(text="Pausa", bind=Pausa)

# Función para cerrar la terminal
def cerrar_terminal():
    os._exit(0)

# Crear un botón que llama a la función cerrar_terminal
button(text="Cerrar programa", bind=cerrar_terminal)

# Función para ajustar el rate según el tamaño de la órbita
def calcular_rate(tamano_orbita):
    base_rate = 30  # Base del rate para órbitas pequeñas
    factor_rate = 500  # Factor de escala para órbitas grandes
    max_rate = 300  # Limitar el rate máximo para evitar problemas de rendimiento
    
    rate_valor = int(base_rate + factor_rate * tamano_orbita)
    rate_valor = min(rate_valor, max_rate)
    
    return rate_valor

rate_valor = calcular_rate(tamano_orbita)


#Utilizamos el método Euler-Cromer para armar la animación

fase = 0
dfase = 0.01
x = 0
dx = 0.01
y = 0
dy = 0.01
dt = 1e3 #paso
t=0
while True:
    rate(rate_valor)
    if corriendo:
        r = est2.pos-est1.pos
        rmag = mag(r)
        rhat = r/rmag
        F = G * est1.masa * est2.masa * r.hat / mag(r)**2
        est1.p = est1.p + F*dt
        est2.p = est2.p - F*dt
        est1.pos = est1.pos + (est1.p/est1.masa) * dt
        est2.pos = est2.pos + (est2.p/est2.masa) * dt
        label1.pos = est1.pos + (est1.p/est1.masa) * dt
        label2.pos = est2.pos + (est2.p/est2.masa) * dt
        vr1_curva.plot(fase,x)
        x=eq_Vr1(eq_tang(kepler.solve(2*np.pi*fase,e)))+dx
        vr2_curva.plot(fase,y)
        y=eq_Vr2(eq_tang(kepler.solve(2*np.pi*fase,e)))+dy
        fase=fase+dfase
        t=t+dt
   
