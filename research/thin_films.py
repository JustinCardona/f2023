import numpy as np
import cmath
import matplotlib.pyplot as plt
from matplotlib import cm


class Medium:
    def __init__(self, mu:complex, epsilon:complex, thickness:float, name='medium'):
        self.m = mu
        self.e = epsilon
        self.d = thickness
        self.name = name


class Wave:
    def __init__(self, a:complex, k:np.ndarray[complex], p='e'):
        self.a = a
        self.k = k
        self.p = p

    def eval(self, m:Medium, t=0, x=np.array([0, 0])):
        phase = self.k@x - cmath.sqrt(self.k@self.k * m.m * m.e) * t
        return self.a * np.exp(1j * phase)
        

def disperse(w:Wave, m_1:Medium, m_2:Medium):
    g = m_2.m * m_2.e / (m_1.m * m_1.e)
    k_T = np.array([w.k[0], cmath.sqrt((w.k@w.k * g - w.k[0]**2))])
    return Wave(w.a, k_T, w.p)


def a_coeff(k_I:np.ndarray[complex], k_T:np.ndarray[float]):
    return (k_T[1] / cmath.sqrt(k_T@k_T)) / (k_I[1] / cmath.sqrt(k_I@k_I))


def b_coeff(m_1:Medium, m_2:Medium):
    return cmath.sqrt((m_1.m * m_2.e) / (m_2.m * m_1.e))


def r(w:Wave, m_1:Medium, m_2:Medium):
    k_T = disperse(w, m_1, m_2).k
    a = a_coeff(w.k, k_T)
    b = b_coeff(m_1, m_2)
    if w.p == 'e':
        return (1 - a*b) / (1 + a*b) # parallel polarization
    else:
        return (a - b) / (a + b) # perpedicular polarization


def t(w:Wave, m_1:Medium, m_2:Medium):
    k_T = disperse(w, m_1, m_2).k
    a = a_coeff(w.k, k_T)
    b = b_coeff(m_1, m_2)
    if w.p == 'e':
        return 2 / (1 + a*b) # parallel polarization 
    else:
        return 2 / (a + b) # perpedicular polarization


def R(w:Wave, M:list[Medium]):
    if len(M)==2:
        return r(w, *M)
    else:
        w_new = disperse(w,  M[0], M[1])
        phase = np.exp(2j*w_new.k[1] * M[1].d)
        r_old = r(w, M[0], M[1])
        r_new = R(w_new, M[1:])
        return (r_old + r_new * phase) / (1 + r_old * r_new * phase)


def T(w:Wave, M:list[Medium]):
    if len(M)==2:
        return t(w, *M)
    else:
        w_new = disperse(w,  M[0], M[1])
        phase = np.exp(1j*w_new.k[1] * M[1].d)
        r_old = r(w, M[0], M[1])
        t_old = t(w, M[0], M[1])
        r_new = R(w_new, M[1:])
        t_new = T(w_new, M[1:])
        return (t_old * t_new * phase) / (1 + r_old * r_new * phase**2)


def plot_coefficients(M:list[Medium], n=1000, K = 1):
    angles = np.linspace(0.001, np.pi / 2, n)
    kx = K * np.append(np.sinh(-np.flip(angles)), np.sin(angles))
    ky = K * np.append(np.cosh(-np.flip(angles)), np.cos(angles))
    re = np.linspace(0, 1, 2*n)
    te = np.linspace(0, 1, 2*n)
    rm = np.linspace(0, 1, 2*n)
    tm = np.linspace(0, 1, 2*n)
    for i in range(len(kx)):
        k = K * np.array([kx[i], ky[i]])
        we = Wave(1, k, 'e')
        wm = Wave(1, k, 'm')
        RE, TE = R(we, M), T(we, M)
        RM, TM = R(wm, M), T(wm, M)
        k_T = disperse(we, M[0], M[-1]).k
        k_I = we.k
        a = a_coeff(k_I, k_T)
        b = b_coeff(M[0], M[-1])
        re[i] = abs(RE)**2
        te[i] = abs(a*b) * abs(TE)**2
        rm[i] = abs(RM)**2
        tm[i] = abs(a*b) * abs(TM)**2

    plt.plot(kx, re, label=r"$R_e$")
    plt.plot(kx, te, label=r"$T_e$")
    plt.plot(kx, rm, label=r"$R_m$")
    plt.plot(kx, tm, label=r"$T_m$")
    plt.xlabel(r"$k_x$")
    plt.legend()
    plt.show()
    plt.close()


def plot_coefficients_3D(M:list[Medium], n=100, K = 1):
    angles = np.linspace(0.001, np.pi / 2, n)
    k_re = np.sin(angles)
    k_y = np.cos(angles)
    k_im = np.sinh(angles - np.pi/4) / np.max(np.sinh(angles - np.pi/4))
    X, Y = np.meshgrid(k_re, k_im)
    re = np.zeros((n, n))
    te = np.zeros((n, n))
    rm = np.zeros((n, n))
    tm = np.zeros((n, n))
    for ix in range(n):
        for iy in range(n):
            k = K * np.array([X[ix][iy] + Y[ix][iy]*1j, k_y[iy]])
            we = Wave(1, K * k, 'e')
            wm = Wave(1, K * k, 'm')
            RE, TE = R(we, M), T(we, M)
            RM, TM = R(wm, M), T(wm, M)
            k_T = disperse(we, M[0], M[-1]).k
            k_I = we.k
            a = a_coeff(k_I, k_T)
            b = b_coeff(M[0], M[-1])
            re[ix][iy] = abs(RE)**2
            te[ix][iy] = abs(a*b) * abs(TE)**2
            rm[ix][iy] = abs(RM)**2
            tm[ix][iy] = abs(a*b) * abs(TM)**2
        
    fig, ax = plt.subplots(2, 2, subplot_kw={"projection": "3d"})    
    color = cm.jet
    surf = ax[0, 0].plot_surface(X, Y , re, cmap=color, linewidth=0, antialiased=False)
    surf = ax[0, 1].plot_surface(X, Y , te, cmap=color, linewidth=0, antialiased=False)
    surf = ax[1, 0].plot_surface(X, Y , rm, cmap=color, linewidth=0, antialiased=False)
    surf = ax[1, 1].plot_surface(X, Y , tm, cmap=color, linewidth=0, antialiased=False)
    ax[0, 0].set_zlabel(r"$R_e$")
    ax[0, 1].set_zlabel(r"$T_e$")
    ax[1, 0].set_zlabel(r"$R_m$")
    ax[1, 1].set_zlabel(r"$T_m$")    
    fig.colorbar(surf)
    for i in [[x, y] for x in [0, 1] for y in [0, 1]]:
        ax[*i].set_xlabel(r"Re($k_x$)")
        ax[*i].set_ylabel(r"Im($k_x$)")
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.show()
    plt.close()


def plot_mult_refl(M:list[Medium], n=1000, d=1e-7):
    ks = np.linspace(1 / (d * np.e**3), np.e**2 / d, n)
    re = np.linspace(0, 1, n)
    te = np.linspace(0, 1, n)
    rm = np.linspace(0, 1, n)
    tm = np.linspace(0, 1, n)
    
    for i in range(len(ks)):
        k = np.array([0, ks[i]])
        we = Wave(1, k, 'e')
        wm = Wave(1, k, 'm')
        RE, TE = R(we, M), T(we, M)
        RM, TM = R(wm, M), T(wm, M)
        k_T = disperse(we, M[0], M[-1]).k
        k_I = we.k
        a = a_coeff(k_I, k_T)
        b = b_coeff(M[0], M[-1])
        re[i] = abs(RE)**2
        te[i] = abs(a*b) * abs(TE)**2
        rm[i] = abs(RM)**2
        tm[i] = abs(a*b) * abs(TM)**2
    plt.plot(ks, re, label=r"$R_e$")
    plt.plot(ks, te, label=r"$T_e$")
    plt.plot(ks, rm, label=r"$R_m$")
    plt.plot(ks, tm, label=r"$T_m$")
    plt.xlabel("$k_y$ (m)")
    plt.legend()
    plt.show()
    plt.close()


def plot_fields(w:Wave, M:list[Medium], n=100):
    rs = [R(w, M[m:]) for m in range(len(M)-1)]
    ts = [T(w, M[:m]) for m in range(2, len(M)+1)]
    ks = [disperse(w, M[0], M[idx]).k for idx in range(1, len(M))]
    A = [ts[idx] * rs[idx+1] * np.exp(2j * ks[idx][1] * M[idx+1].d) / (1 - rs[idx] * rs[idx+1] * np.exp(2j * ks[idx][1] * M[idx+1].d)) for idx in range(len(ts)-1)]

    ws = [Wave(w.a*rs[0], w.k, w.p)] + [Wave(w.a * A[idx], ks[idx], w.p) for idx in range(len(A))] + [Wave(w.a*(ts[-1]), ks[-1], w.p)]
    ds = np.cumsum([0] + [m.d for m in M])
    max_I = 0
    for idx in range(len(ws)):
        domain = np.linspace(ds[idx], ds[idx+1], n)
        field = list(map(lambda y: 0.5 * np.abs(M[idx].m * M[idx].e * (ws[idx].eval(M[idx], 0, [0, y-ds[idx]]))**2), domain))
        max_I = max(max_I,  max(field))
        plt.plot(domain, field, label=M[idx].name)
    plt.vlines(ds[1:-1], 0, max_I, colors='black')
    plt.legend()
    plt.show()
    plt.close()
    
# TESTING
vaccuum = Medium(1, 1, 1e-7, 'vaccuum')
water = Medium(1, 1.33**2, 3e-9, 'water')
glass = Medium(1, 2.2, 1e-9, 'glass')
dielectric = Medium(1, 10, 1e-7, 'dielectric')
smtg = Medium(1, 1.82+0.08j, 1e-7, 'smtg')
silver = Medium(1, -18.2+0.5j, 5.3e-8, 'silver')
gold = Medium(1, -11.6+1.2j, 5.0e-8, 'gold')
# m_test = [vaccuum, glass, smtg, vaccuum]
m_test = [vaccuum, dielectric, silver, dielectric]

k_test = (1e7) * np.array([np.sqrt(2) / 2, np.sqrt(2) / 2])
we_test = Wave(1, k_test, 'e')
wm_test = Wave(1, k_test, 'm')

plot_coefficients(m_test)
plot_coefficients_3D(m_test)
plot_mult_refl(m_test)
plot_fields(we_test, m_test)