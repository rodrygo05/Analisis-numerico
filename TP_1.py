import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import metodos_numericos as met
import funciones as fn

##arreglo de funciones 
funciones=[fn.funcion1 ,fn.funcion2, fn.funcion3]
derivada_pri=[fn.funcion1p ,fn.funcion2p, fn.funcion3p]
derivada_seg=[fn.funcion1pp ,fn.funcion2pp, fn.funcion3pp]

titulos=["funcion 1","funcion 2","funcion 3"] #arreglo de cadenas para cambiar titulos 

# Paramentros de la configuracion
tolerancia = 1e-13
nMax = 50
parametros_a  = [0, 0, 0]
parametros_b  = [2, 2, 2]
parametros_p0 = [1.250000, 1.125000, 1.250000]
parametros_p1 = [1.375000, 1.187500, 1.375000]

for i in range(3):
    a = parametros_a[i]
    b = parametros_b[i]
    p0 = parametros_p0[i]
    p1 = parametros_p1[i]

    x = np.linspace(0,2, 1000)
    Yfuncion = funciones[i](x)
    plt.figure()
    plt.plot(x, Yfuncion, 'k', lw=1)    
    nombreGrafico1 = titulos[i]
    plt.title(nombreGrafico1)
    plt.grid(True)

    # BISECCION
    print("Resultados biseccion " +str(titulos[i]))
    raiz_biseccion,nIteraciones_biseccion,historiaRaices_biseccion = met.biseccion(funciones[i],a,b,tolerancia,nMax)
    # SECANTE
    print("Resultados secante " +str(titulos[i]))
    raiz_sec,nIteraciones_sec,historiaRaices_sec = met.secante(funciones[i],p0,p1,tolerancia,nMax)
    # NEWTON RAPHSON
    print("Resultados Newton Raphson " +str(titulos[i]))
    raiz_NR,nIteraciones_NR,historiaRaices_NR = met.NR(funciones[i],derivada_pri[i],p0,tolerancia,nMax)
    # NEWTON RAPHSON MODIFICADO
    print("Resultados Newton Raphson modificado " +str(titulos[i]))
    raiz_NR_mod,nIteraciones_NR_mod,historiaRaices_NR_mod = met.NR_mod(funciones[i],derivada_pri[i],derivada_seg[i],p1,tolerancia,nMax)  
   
    plt.figure()
    plt.plot(historiaRaices_biseccion[:,0],historiaRaices_biseccion[:,1],'-',lw=2,label='Biseccion')
    plt.plot(historiaRaices_sec[:,0], historiaRaices_sec[:,1],'-',lw=2, label='Secante')
    plt.plot(historiaRaices_NR[:,0], historiaRaices_NR[:,1],'-',lw=2, label='Newton Raphson ')
    plt.plot(historiaRaices_NR_mod[:,0],historiaRaices_NR_mod[:,1],'-',lw=2,label='Newton Raphson  modificado')
    
    plt.xlabel('Paso [n]')
    plt.title('Raiz estimada ' +str(titulos[i]))
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()
       
    ordenConvergenciaBiseccion = met.estimarOrdenConvergencia(historiaRaices_biseccion,nIteraciones_biseccion)
    
    ordenConvergencia_sec = met.estimarOrdenConvergencia(historiaRaices_sec,nIteraciones_sec)
    
    ordenConvergencia_NR = met.estimarOrdenConvergencia(historiaRaices_NR,nIteraciones_NR)
    
    ordenConvergencia_NR_mod = met.estimarOrdenConvergencia(historiaRaices_NR_mod,nIteraciones_NR_mod)
    
    plt.figure()
    plt.plot(ordenConvergenciaBiseccion[:,0],ordenConvergenciaBiseccion[:,1],'-',lw=2,label='Biseccion')
    plt.plot(ordenConvergencia_sec[:,0], ordenConvergencia_sec[:,1],'-', lw=2, label='secante')
    plt.plot(ordenConvergencia_NR[:,0], ordenConvergencia_NR[:,1],'-', lw=2, label='newton raphson')
    plt.plot(ordenConvergencia_NR_mod[:,0],ordenConvergencia_NR_mod[:,1],'-',lw=2,label='newton raphson modificado')
    
    plt.xlabel('Paso [n]')
    plt.ylabel('alfa')
    plt.title('Orden de convergencia ' +str(titulos[i]))
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()
    
    Lambda_biseccion=met.contante_asintotica(ordenConvergenciaBiseccion,historiaRaices_biseccion,nIteraciones_biseccion)
    
    Lambda_sec=met.contante_asintotica(ordenConvergencia_sec,historiaRaices_sec,nIteraciones_sec)
        
    Lambda_NR=met.contante_asintotica(ordenConvergencia_NR,historiaRaices_NR,nIteraciones_NR)
        
    Lambda_NR_mod=met.contante_asintotica(ordenConvergencia_NR_mod,historiaRaices_NR_mod,nIteraciones_NR_mod)
        
    plt.figure()
    plt.plot(Lambda_biseccion[:,0],Lambda_biseccion[:,1],'-',lw=2,label='Biseccion')
    plt.plot(Lambda_NR[:,0], Lambda_NR[:,1],'-', lw=2, label='secante')
    plt.plot(Lambda_sec[:,0], Lambda_sec[:,1],'-', lw=2, label='newton raphson')
    plt.plot(Lambda_NR_mod[:,0],Lambda_NR_mod[:,1],'-',lw=2,label='newton raphson modificado')
        
    plt.xlabel('Paso [n]')
    plt.ylabel('lambda')
    plt.title('Constante asintotica ' +str(titulos[i]))
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()
    
print("VALORES CALCULADOS EMPLEANDO LA BIBLIOTECA SCIPY\n")

print("Bisección:")
for i in range(3):
    a = parametros_a[i]
    b = parametros_b[i]
    print("\tLa raíz de la %s es % .15f\n"%(titulos[i],optimize.bisect(funciones[i],a,b)))

print("Newton-Raphson Ordinario / Secante:")
for i in range(3):
    p0 = parametros_p0[i]
    try:
        print("\tLa raíz de la %s es % .15f\n"%(titulos[i],optimize.newton(funciones[i],p0)))
    except RuntimeError as e:
     print("\tAl calcular la raiz de la %s:\n\t"%(titulos[i]) + str(e) + "\n")

print("Newton-Raphson Modificado:")
for i in range(3):
    p0 = parametros_p0[i]
    try:
        print("\tLa raíz de la %s es % .15f\n"%(titulos[i],optimize.newton(funciones[i],p0,derivada_pri[i])))
    except RuntimeError as e:
     print("\tAl calcular la raiz de la %s:\n\t"%(titulos[i]) + str(e) + "\n")

