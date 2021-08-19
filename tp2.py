import metodos as mt
import matplotlib.pyplot as plt
import numpy as np

# Funciones del sistema de ecuaciones

def f(t, u1, u2):
    return u2

def g(t, u1, u2):
    return -((b/m)*u2 + (gravedad/l)*u1)
    
# Funciones auxiliares

def energia(historia):
    largo=len(historia)
    E=np.zeros((largo,1))
    for i in range(0,largo):
        E[i]=m*gravedad*l*(1 - np.cos(historia[i][1]))+0.5*m*(l*historia[i][2])**2
    return E
    
def imprimir(historia):
    largo=len(historia)
    if(largo>10):
        hist_prim=historia[0:5]
        hist_ult=historia[-5:]
        print("------------------------------------------")
        print("|      t     |      u1     |     u2      |")
        print("------------------------------------------")
        for punto in hist_prim:
            print("|  % .6f | % .3e  | % .3e  |"%(punto[0], punto[1], punto[2]))
        print("|     .      |      .      |      .      |")
        print("|     .      |      .      |      .      |")
        print("|     .      |      .      |      .      |")
        for punto in hist_ult:
            print("| % .6f | % .3e  | % .3e  |"%(punto[0], punto[1], punto[2]))
        print("------------------------------------------")
        print("\n")
    elif(largo<=10):
        print("---------------------------'--------------")
        print("|     t      |      u1     |     u2      |")
        print("---------------------------'--------------")
        for punto in historia:
            print("| % .6f  | % .3e  | % .3e  |"%(punto[0], punto[1], punto[2]))
        print("------------------------------------------")
        
# Parametros por defecto
        
gravedad = 9.81
m = 1
l = 1
b = 0 # 0.5
theta0 = 0.5235987756
thetax0 = 0 #1.745329252
tiempoi = 0
tiempof = 20
paso = 0.2

# Entrada de usuario

print("CALCULADORA PENDULO AMORTIGUADO ")
print("¿Usar valores por defecto? (s/n)")
if(input() != "s"):
    print("Ingresar masa [Kg]: ", end="")
    m = float(input())
    print("Ingresar longitud del hilo [m]: ", end="")
    l = float(input())
    print("Ingresar coeficiente de amortiguamiento [N.s/m]: ", end="")
    b = float(input())
    print("Ingresar angulo inicial [rad]: ", end="")
    theta0 =  float(input())
    print("Ingresar la velocidad angular inicial [rad/s]: ", end="")
    theta0x =  float(input())
    print("Ingresar tiempo inicial [s]: ", end="")
    tiempoi = float(input())
    print("Ingresar tiempo final [s]: ", end="")
    tiempof = float(input())
    print("Ingresar paso [s]: ", end="")
    paso = float(input())
print("\n")
    
historia_e = mt.resolver_ec_dif_segundo_orden(mt.euler, f, g, paso,\
                                            tiempoi, theta0, thetax0, tiempof)
    
historia_rk4 = mt.resolver_ec_dif_segundo_orden(mt.runge_kutta_4, f, g, paso,\
                                            tiempoi, theta0, thetax0, tiempof)
# Salida a consola
    
print("EULER")
imprimir(historia_e)
 
print("RK4")
imprimir(historia_rk4)

historia_e = np.array(historia_e)
historia_rk4 = np.array(historia_rk4)

E_e=energia(historia_e)
E_rk4=energia(historia_rk4)

# Gráfico

fig, axs = plt.subplots(3,2,sharex=True) 

axs[0][0].plot(historia_e[:,0], historia_e[:,1],"b-")
axs[0][0].set_title('RK1')
axs[0][0].set_ylabel('Pos[rad]')
axs[0][0].grid(True)

axs[0][1].plot( historia_rk4[:,0], historia_rk4[:,1],"r-")
axs[0][1].set_title('RK4')
axs[0][1].grid(True)

axs[1][0].plot(historia_e[:,0], historia_e[:,2],"b-")
axs[1][0].set_ylabel('Vel[rad/s]')
axs[1][0].grid(True)

axs[1][1].plot( historia_rk4[:,0], historia_rk4[:,2],"r-")
axs[1][1].grid(True)

axs[2][0].plot(historia_e[:,0], E_e,"b-")
axs[2][0].grid(True)
axs[2][0].set_ylabel('E[J]')
axs[2][0].set_xlabel('Tiempo[s]')

axs[2][1].plot( historia_rk4[:,0], E_rk4,"r-")
axs[2][1].grid(True)
axs[2][1].set_xlabel('Tiempo[s]')

fig.tight_layout()
fig.savefig("fig1.svg" ,orientation='portrait')
plt.show()
