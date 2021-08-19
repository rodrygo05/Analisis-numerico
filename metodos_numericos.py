# -*- coding: utf-8 -*-
import numpy as np

def secante(f,semilla0,semilla1,TOL,N0):
    historia = np.zeros((N0,2))
    p_n1=semilla1  #Pn-1
    p_n2=semilla0 #Pn-2
    print("------------------------------------------------------------------------")
    print("|  n   |         p(n)        |   p_(n-1)   |   p_(n-2)   |      e      |")
    print("------------------------------------------------------------------------")
    print("| %4d |         -           | % f   | % f   |      -      |"%(0,p_n1,p_n2))
    for i in range(1,N0):
        p_n= p_n1- f(p_n1) * (p_n1 - p_n2) / ( f(p_n1)-f(p_n2) )
        historia[i] = (i,p_n)
        error = p_n-p_n1
        if np.abs(error) <= TOL:
            historia = historia[:i+1]
            if(error==0):
                print("| %4d | % .15f  | % f   | % f   | % .3e |"%(i,p_n,p_n1,p_n2,2.225e-308))
            else:
                print("| %4d | % .15f  | % f   | % f   | % .4e |"%(i,p_n,p_n1,p_n2,error))
            print("------------------------------------------------------------------------")
            print('\n')
            return p_n, i,historia                  
        else:
            print("| %4d | % .15f  | % f   | % f   | % .4e |"%(i,p_n,p_n1,p_n2,error))
            p_n2=p_n1 
            p_n1=p_n
            
    
    print("Luego de %d iteraciones no se obtuvo un resultado." %N0)
    historia = historia[:i+1]
    return p_n, i,historia  
    


        
def NR_mod(f,fp,fpp,semilla,TOL,N0):
    print("-------------------------------------------")
    print("|  n   |         p          |      e      |")
    print("-------------------------------------------")
    print("| %4d | % .15f |      -      |"%(0,semilla))
   
    historia = np.zeros((N0,2))
    p_anterior=semilla
    historia[0] = (0,p_anterior)
    for i in range(1,N0):
        
        p_actual=p_anterior- (f(p_anterior)*fp(p_anterior)) / ((fp(p_anterior)**2)-(f(p_anterior)*fpp(p_anterior)))
        error = (p_actual-p_anterior);
        historia[i] = (i,p_actual)
        if(error==0):
            print("| %4d | % .15f | % .3e |"%(i,p_actual,2.225e-308))
        else:    
            print("| %4d | % .15f | % .4e |"%(i,p_actual,error))
        if abs(error) <= TOL and i > 1:
            historia = historia[:i+1]
            print("-------------------------------------------")
            print('\n')
            return p_actual, i,historia                  
        
        p_anterior=p_actual 
    historia = historia[:i+1]
    return p_actual, i,historia  


def NR(f,fp,semilla,TOL,N0):
    print("-------------------------------------------")
    print("|  n   |         p          |      e      |")
    print("-------------------------------------------")
    print("| %4d | % .15f |      -      |"%(0,semilla))
   
    historia = np.zeros((N0,2))
    p_anterior=semilla
    historia[0] = (0,p_anterior)
    for i in range(1,N0):
        
        p_actual=p_anterior-(f(p_anterior)/fp(p_anterior))
        error = (p_actual-p_anterior)
        historia[i] = (i,p_actual)
        if(error==0):
            print("| %4d | % .15f | % .3e |"%(i,p_actual,2.225e-308))
        else:    
            print("| %4d | % .15f | % .4e |"%(i,p_actual,error))
        if abs(error) < TOL and i > 1:
            historia = historia[:i+1]
            print("-------------------------------------------")
            print('\n')
            return p_actual, i,historia                  
        
        p_anterior=p_actual 
    historia = historia[:i+1]
    return p_actual, i,historia

           
def biseccion(f,a,b,TOL,N0):

    if f(a)*f(b) >= 0:
        print("No puede existir una raiz en el intervalo dado.")
        return None

    historia = np.zeros((N0,2))
    p_anterior = a
    historia[0] = (0,p_anterior)
    
    print("-----------------------------------------------------------------------------------")
    print("|  n   |     a     |      b    |          p         |      f(p)     |      e      |")
    print("-----------------------------------------------------------------------------------")
    print("| %4d | % f | % f |         -          |       -       |       -     |"%(0,a,b))
    for i in range(1,N0):
        p = (a+b)/2
        historia[i] = (i,p)
        fun=f(p)
        error=abs(fun-f(p_anterior)) 
        if error <= TOL and i > 1:
            historia = historia[:i+1]
            print("| %4d | % f | % f | % .15f | % e | % .4e |"%(i,a,b,p,fun,error))
            print("-----------------------------------------------------------------------------------")
            print('\n')
            return p, i,historia
        
        if f(a)*f(p) > 0: # No tengo nada a la izquierda.
            a = p
        else: # No tengo nada a la derecha.
            b = p
        print("| %4d | % f | % f | % .15f | % e | % .4e |"%(i,a,b,p,fun,error))  
        p_anterior = p
        
    print("No convergià¸£à¸“.")
    historia = historia[:i+1]
    return p, i,historia

def estimarOrdenConvergencia(historiaRaices,nInteraciones):
        
        alfa = np.zeros((nInteraciones-1,2))
        
        for n in range(3-1,nInteraciones-1):
            e_n_mas_1 = np.abs(historiaRaices[n+1][1]-historiaRaices[n][1])
            e_n = np.abs(historiaRaices[n][1] - historiaRaices[n-1][1])
            e_n_menos_1 = np.abs(historiaRaices[n-1][1]-historiaRaices[n-2][1])
            a=np.log( np.abs(e_n_mas_1 / e_n ))
            b=np.log( np.abs(e_n / e_n_menos_1 ))
            alfa[n] = n, a/b
            
        return alfa
        
def contante_asintotica(alfa,historiaRaices,nIteraciones):
    
        Lambda = np.zeros((nIteraciones-1,2))
        
        for n in range(2-1,nIteraciones-1):
            e_n_mas_1 = np.abs(historiaRaices[n+1][1]-historiaRaices[n][1])
            e_n = np.abs(historiaRaices[n][1] - historiaRaices[n-1][1])
            Lambda[n] = n, e_n_mas_1/(e_n**(alfa[n][1]))
            
        return Lambda
