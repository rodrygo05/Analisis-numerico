# Metodos de resolucion de ecuaciones diferenciales
# Cada uno retorna un avance del método homónimo.

def euler(f, g, paso, t, u1, u2):
    return (paso*f(t, u1, u2), paso*g(t, u1, u2))

def runge_kutta_4(f, g, paso, t, u1, u2):
    k1 = f(t, u1, u2)
    l1 = g(t, u1, u2)
    
    k2 = f(t + paso/2, u1 + (paso/2)*k1, u2 + (paso/2)*l1)
    l2 = g(t + paso/2, u1 + (paso/2)*k1, u2 + (paso/2)*l1)
    
    k3 = f(t + paso/2, u1 + (paso/2)*k2, u2 + (paso/2)*l2)
    l3 = g(t + paso/2, u1 + (paso/2)*k2, u2 + (paso/2)*l2)
    
    k4 = f(t + paso, u1 + paso*k3, u2 + paso*l3)
    l4 = g(t + paso, u1 + paso*k3, u2 + paso*l3)
    
    u1x = (paso/6)*(k1 + 2*k2 + 2*k3 + k4)
    u2x = (paso/6)*(l1 + 2*l2 + 2*l3 + l4)
    
    return (u1x, u2x)


# Resuelve ecuaciones de la forma:
#
# du1 = f(u1, u2, t)
# dt
#
# du2 = g(u1, u2, t)
# dt
#
# Con u1 y u2 conocidos en un t0 
    
def resolver_ec_dif_segundo_orden\
(metodo, f, g, paso, t0, u10, u20, tf):
    historia = []
    ti = t0
    u1i = u10
    u2i = u20

    while(ti <= tf):
        historia.append((ti,u1i,u2i))
        u1x, u2x = metodo(f, g, paso, t0, u1i, u2i)
        u1i += u1x
        u2i += u2x
        ti += paso
    return historia