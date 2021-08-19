from scipy import optimize
import funciones as fn

funciones=[fn.funcion1 ,fn.funcion2, fn.funcion3]
derivada_pri=[fn.funcion1p ,fn.funcion2p, fn.funcion3p]
titulos=["funcion 1","funcion 2","funcion 3"]

a = 1
b = 2
nMax = 100
p0 = 1.125000

print("VALORES CALCULADOS EMPLEANDO LA BIBLIOTECA SCIPY\n")

print("Bisección:")
for i in range(3):
    print("\tLa raíz de la %s es % .15f\n"%(titulos[i],optimize.bisect(funciones[i],a,b)))

print("Newton-Raphson Ordinario / Secante:")
for i in range(3):
    try:
        print("\tLa raíz de la %s es % .15f\n"%(titulos[i],optimize.newton(funciones[i],p0)))
    except RuntimeError as e:
     print("\tAl calcular la raiz de la %s:\n\t"%(titulos[i]) + str(e) + "\n")


print("Newton-Raphson Modificado:")
for i in range(3):
    try:
        print("\tLa raíz de la %s es % .15f\n"%(titulos[i],optimize.newton(funciones[i],p0,derivada_pri[i])))
    except RuntimeError as e:
     print("\tAl calcular la raiz de la %s:\n\t"%(titulos[i]) + str(e) + "\n")

              