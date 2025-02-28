##Programme de simulation d'oscillations
#FRINZI Mathis

#Simulation tkinter
try:
    from tkinter import *
except:
    from Tkinter import *

from math import *
import time

def x(a):
    return a[0]

def v(a):
    return a[1]

def a(x):
    return x[2]

class Simu():
    ''' Stocke les variables propre à la simulation : constante de la tour etc...'''
    def __init__(self,fenetre,th):
        self.fenetre=fenetre
        self.L=0
        self.M=0
        self.M1=0
        self.K=0
        self.K1=0
        self.KP=0
        self.H=0
        self.H1=0
        self.H2=0
        self.Fex = lambda x:100
        self.theorie = th
        self.Eul = Euler(self)
        self.EulS = EulerSansPendule(self)
        self.tour = Tour(self.fenetre)
        self.tourS = ToursansPendule(self.fenetre)

class EulerSansPendule():
    '''Méthode d'Euler : construction d'une solution à l'équation différentielle'''
    def __init__(self,sim):
        self.X=[0,0,0]
        self.simu = sim
        if self.simu.theorie==1:
            self.PsiX=self.PsiX1
        if self.simu.theorie==0:
            self.PsiX=self.PsiX0

    def PsiX0(self,t,h):
        ''' Fonction qui à un couple (pos,vit) renvoie le couple(pos',vit')=(vit,acc)'''
        X=self.X
        acX = self.simu.Fex(t)+1/(self.simu.M) * (-self.simu.K*x(X)-self.simu.H*v(X))
        return v(X),acX

    def PsiX1(self,t,h):
        ''' Fonction qui à un couple (pos,vit) renvoie le couple(pos',vit')=(vit,acc)'''
        X=self.X
        dxs = lambda t:(self.simu.Fex(t-h)-self.simu.Fex(t))/h
        ddxs = (dxs(t-h)-dxs(t))/h
        acX = 1/(self.simu.M) * (-self.simu.K*(x(X)-self.simu.Fex(t)) -self.simu.H*(v(X)-dxs(t)))
        return v(X),acX

    def PointSuivant(self,t,pas=0.001):
        '''Fonction qui plot un nouveau point à partir des coordonnées du précédent'''
        ancienX = self.X.copy()
        #NouveauX
        CoupleX = self.PsiX(t,pas)
        accelerationX = CoupleX[1]
        vitesseX = v(ancienX) + pas*CoupleX[1]
        positionX = x(ancienX) + pas*CoupleX[0] + pas**2/2 * CoupleX[1]
        self.X = [positionX,vitesseX,accelerationX]

class Euler():
    '''Méthode d'Euler : construit la solution de l'équation différentielle'''
    def __init__(self,sim):
        self.X=[0.5,0,0]
        self.U=[0,0,0]
        self.simu = sim
        if self.simu.theorie==1:
            self.PsiX=self.PsiX1
            self.PsiU=self.PsiU1
        if self.simu.theorie==0:
            self.PsiX=self.PsiX0
            self.PsiU=self.PsiU0

    def PsiX0(self,t,h):
        ''' Equation différentielle de la tour'''
        X=self.X
        U=self.U
        u = self.U[0]
        M=self.simu.M
        M1=self.simu.M1
        K=self.simu.K
        H=self.simu.H
        L=self.simu.L
        K1=self.simu.K1
        KP=self.simu.KP
        acX = self.simu.Fex(t)+1/(M+M1) * (-K*x(X)-H*v(X)- M1*L * a(U)* cos(u) + M1*L*(v(U))**2 * sin(u)) + (K1+KP)*L*sin(u)
        return v(X),acX

    def PsiX1(self,t,h):
        ''' Fonction qui à un couple (pos,vit) renvoie le couple(pos',vit')=(vit,acc)'''
        X=self.X
        U=self.U
        u = self.U[0]
        self.simu.K1 = self.simu.M1*9.81*(self.simu.L**2+self.simu.Eul.U[0]**2)**0.5/(self.simu.L**2) #formule  de théorie
        M=self.simu.M
        M1=self.simu.M1
        K=self.simu.K
        H=self.simu.H
        H1=self.simu.H1
        L=self.simu.L
        K1=self.simu.K1
        KP=self.simu.KP
        w0T = (K/M)**0.5
        e = (K1+KP)/K
        P = H1/H
        zT = H/(2*((K1+KP)*M)**0.5)
        dxs = (self.simu.Fex(t)-self.simu.Fex(t-h))/h
        acX = (-w0T**2) *(1+e)*x(X) -2*w0T*zT*(1+P)*v(X) + e*(w0T**2)*x(U) + 2*P*w0T*zT*v(U) + (w0T**2)*self.simu.Fex(t)+2*w0T*zT*dxs
        return v(X),acX

    def PsiU0(self,t,h):
        '''Fonction qui à un couple (pos(U),vit(U)) renvoie (vit(U),acc(U))'''
        X = self.X
        U = self.U
        u = U[0]
        if self.simu.M1!=0 and self.simu.L != 0:
            self.simu.K1 = self.simu.M1*9.81*(self.simu.L**2+self.simu.Eul.U[0]**2)**0.5/(self.simu.L**2)
            acU = 1/self.simu.L * (-self.simu.M1*cos(u)*a(X)- (self.simu.K1+self.simu.KP)*self.simu.L*sin(u)*cos(u) - self.simu.M1*9.81*sin(u) - self.simu.H1*v(U) - self.simu.H2*v(U))
        else:
            return 0,0
        return v(U), acU

    def PsiU1(self,t,h):
        '''Fonction qui à un couple (pos(U),vit(U)) renvoie (vit(U),acc(U))'''
        X = self.X
        U = self.U
        u = U[0]
        if self.simu.M1!=0 and self.simu.L != 0:
            w0P = ((self.simu.K1+self.simu.KP)/self.simu.M1)**0.5
            zP = 0.5*self.simu.H1/(((self.simu.K1+self.simu.KP)*self.simu.M1)**0.5)
            acU = + (w0P**2) * (x(X)-x(U)) + 2*zP*w0P*(v(X)-v(U))
        else:
            return 0,0
        return v(U), acU

    def PointSuivant(self,t,pas=0.0001):
        '''Fonction qui plot un nouveau point à partir des coordonnées du précédent'''
        ancienX = self.X.copy()
        ancienU = self.U.copy()
        #NouveauX
        CoupleX = self.PsiX(t,pas)
        accelerationX = CoupleX[1]
        vitesseX = v(ancienX) + pas*CoupleX[1]
        positionX = x(ancienX) + pas*CoupleX[0] + pas**2/2 * CoupleX[1]
        #NouveauU
        CoupleU = self.PsiU(t,pas)
        accelerationU = CoupleU[1]
        vitesseU = v(ancienU) + pas*CoupleU[1]
        positionU = x(ancienU) + pas*CoupleU[0] + pas**2/2 * CoupleU[1]

        self.X,self.U= [positionX,vitesseX,accelerationX],[positionU,vitesseU,accelerationU]

class Tour():
    '''Stocke les variables de la tour, fais la transition entre l'interface graphique et les données réelles'''
    def __init__(self,fenetre):
        self.fenetre = fenetre
        self.x=0
        self.u=0
        self.list=()

    def unpack(self):
        for i in self.list:
            self.fenetre.canvas.delete(i)

    def pack(self):
        self.unpack()
        a=[]
        # 50 pixels au lieu de 12.5cm  rapport * 4
        if self.fenetre.simu.theorie == 1:
            a.append(self.fenetre.canvas.create_rectangle(250-50 + self.x,100,self.x+250+50,200,fill='#E5E8E8'))
            a.append(self.fenetre.canvas.create_rectangle(300+self.u-5,140,self.u+300+5,160,fill='brown'))
            a.append(self.fenetre.canvas.create_line(300+self.u-5,145,self.x+250+50,145))
            a.append(self.fenetre.canvas.create_line(300+self.u-5, 155,self.x+250+50,155))
            a.append(self.fenetre.canvas.create_line(self.x+250-50, 105, self.fenetre.simu.Fex(self.fenetre.t0)+250-75, 105 ))
            a.append(self.fenetre.canvas.create_line(self.x+250-50, 195, self.fenetre.simu.Fex(self.fenetre.t0)+250-75, 195 ))
            a.append(self.fenetre.canvas.create_line( self.fenetre.simu.Fex(self.fenetre.t0)+250-75,100,  self.fenetre.simu.Fex(self.fenetre.t0)+250-75,200))
        if self.fenetre.simu.theorie == 0:
            a.append(self.fenetre.canvas.create_polygon(125-50 + self.x,100,125-50,500,125+50,500,self.x+125+50,100,fill='#E5E8E8'))
            a.append(self.fenetre.canvas.create_line(125-50 + self.x,100,125-50,500))
            a.append(self.fenetre.canvas.create_line(125+50,500,self.x+125+50,100))
            a.append(self.fenetre.canvas.create_line(125-50 + self.x,100,self.x+125+50,100))
            a.append(self.fenetre.canvas.create_line(125+self.x,100, 125+400*self.fenetre.simu.L*sin(self.u/400),100+400*self.fenetre.simu.L*cos(self.u/400), fill='red'))
            a.append(self.fenetre.canvas.create_line(0,490-1 *self.fenetre.simu.Fex(self.fenetre.t0)/10, 500, 490+1*self.fenetre.simu.Fex(self.fenetre.t0)/10, fill = 'brown'))
        self.list = a

class ToursansPendule():
    '''Stocke les variables de la tour, fais la transition entre l'interface graphique et les données réelles'''
    def __init__(self,fenetre):
        self.fenetre = fenetre
        self.x=0
        self.u=0
        self.list=()
    def unpack(self):
        for i in self.list:
            self.fenetre.canvas.delete(i)
    def pack(self):
        self.unpack()
        a=[]
        # 50 pixels au lieu de 12.5cm  rapport * 4
        if self.fenetre.simu.theorie == 1:
            a.append(self.fenetre.canvas.create_rectangle(250-50 + self.x,400,self.x+250+50,500,fill='#E5E8E8'))
            a.append(self.fenetre.canvas.create_line(self.x+250-50, 405, self.fenetre.simu.Fex(self.fenetre.t0)+250-75, 405 ))
            a.append(self.fenetre.canvas.create_line(self.x+250-50, 495, self.fenetre.simu.Fex(self.fenetre.t0)+250-75, 495 ))
            a.append(self.fenetre.canvas.create_line( self.fenetre.simu.Fex(self.fenetre.t0)+250-75,400,  self.fenetre.simu.Fex(self.fenetre.t0)+250-75,500))
        if self.fenetre.simu.theorie == 0:
            a.append(self.fenetre.canvas.create_polygon(375-50 + self.x,100,375-50,500,375+50,500,self.x+375+50,100,fill='#E5E8E8'))
            a.append(self.fenetre.canvas.create_line(375-50 + self.x,100,375-50,500))
            a.append(self.fenetre.canvas.create_line(375+50,500,self.x+375+50,100))
            a.append(self.fenetre.canvas.create_line(375-50 + self.x,100,self.x+375+50,100))
        self.list = a

class Bode():
    def __init__(self,fenetre):
        ''' Génére l'interface diagramme de Bode '''
        self.fenetreBode = Tk()
        self.fenetre=fenetre
        self.fenetreBode.title('Diagramme de Bode issu de la théorie 1')
        self.panedBode=PanedWindow(self.fenetreBode)
        self.simu = self.fenetre.simu
        if self.simu.theorie==1:
            self.fctDeTrsfrt = self.fctDeTrsfrt1
            self.fctSanspendule = self.fctSanspendule1
        if self.simu.theorie==0:
            self.fctDeTrsfrt = self.fctDeTrsfrt0
            self.fctSanspendule = self.fctSanspendule0
        #Initialisation de l'interface diagramme de Bode
        self.canvas = Canvas(self.panedBode, height=700, width=600,bg='white')
        self.canvas.create_text(550,410, text = "tour avec pendule", fill = 'blue')
        self.canvas.create_text(550,430, text = "tour sans pendule", fill = 'purple')
        self.canvas.create_text(530,450, text = "pendule (référentiel tour)", fill = 'red')
        self.panedBode.add(self.canvas)
        self.panedMenu= PanedWindow(self.panedBode, orient=VERTICAL)
        self.rapp=StringVar(self.fenetreBode,'10')
        self.rapE = Entry(self.panedMenu, text=self.rapp)
        self.panedMenu.add(self.rapE)
        self.bouton = Button(self.panedMenu,text = "Modifier l'échelle (w en rad/s)", command=self.dessiner)
        self.panedMenu.add(self.bouton)
        self.bouton2 = Button(self.panedMenu,text = "Paramètres idéaux", command=self.modific)
        self.panedMenu.add(self.bouton2)
        self.canvasParametres = Canvas(self.panedMenu, width=150,height=400, bg='white')
        self.panedMenu.add(self.canvasParametres)
        self.panedBode.add(self.panedMenu)
        self.ideal = False
        self.panedBode.pack()
        self.w=[]
        self.f=[]
        self.f2=[]
        self.liste=[]
        self.listeP=[]
        self.listespe=[]
        self.dessiner()
        self.after()
        self.fenetreBode.mainloop()

    def fctDeTrsfrt0(self,w):
        '''Diagramme de Bode de la tour avec pendule dans la théorie 0'''
        w0P = (9.81/self.simu.L)**0.5
        w0T = (self.simu.K/self.simu.M)
        Q = self.simu.M/self.simu.H * w0T
        wKP = ((self.simu.K1+self.simu.KP)/self.simu.M1)**0.5
        p = self.simu.M1/(self.simu.M1+self.simu.M)
        return w**2 * ((1-p)*w**2-(1+p)*wKP**2-w0P**2 - 1j*w*w0P/Q)/( (1-p)*(w0T**2-w**2 + 1j*w*w0T/Q)*(w**2-wKP**2- w0P**2-1j*w*w0P/Q) + w**2*p*(w**2+wKP**2))

    def fctDeTrsfrt1(self,w):
        '''Diagramme de Bode de la tour avec pendule dans la théorie 1'''
        w0P = ((self.simu.K1+self.simu.KP)/self.simu.M1)**0.5
        w0T = (self.simu.K/self.simu.M)**0.5
        e =(self.simu.K1+self.simu.KP)/self.simu.K
        P = self.simu.H1/self.simu.H
        zP = self.simu.H1/(2*((self.simu.K1+self.simu.KP)*self.simu.M1)**0.5)
        zT = self.simu.H/(2*(self.simu.K*self.simu.M)**0.5)
        return ((1j*w/w0P)**2+2*zP*(1j*w)/w0P + 1)*(2*zT*1j*w/w0T+1)/(((1j*w/w0P)**2 + 2*zP*1j*w/w0P + 1)*( (1j*w/w0T)**2 + 2*zT*1j*w/w0T + 1) + ((1j*w/w0P)**2)*(2*zT*P*1j*w/w0T + e))

    def fctSanspendule1(self,w):
        '''Diagramme de Bode de la tour sans pendule dans la théorie 1'''
        w0T = (self.simu.K/self.simu.M)**0.5
        zT = self.simu.H/(2*(self.simu.K*self.simu.M)**0.5)
        return (2*zT*1j*w/w0T + 1)/((1j*w/w0T)**2 + 2*zT*1j*w/w0T + 1)
        # return (w0T**2 + 2*w0T*zT*1j*w)/(w0T**2 - w**2 + 2*w0T*zT*1j*w )
        # return (1+2j*w*zT/w0T)/((w/w0T)**2 + 2*zT*1j*w/w0T + 1)

    def fctSanspendule0(self,w):
        '''Diagramme de Bode de la tour sans pendule dans la théorie 0'''
        w0T = (self.simu.K/self.simu.M)**0.5
        QT = self.simu.M1/self.simu.H1 *w0T
        return (w**2)/(-w**2+w0T**2 +1j*w0T*w/QT)

    def fctPendule(self,w):
        '''Diagramme de Bode du pendule'''
        w0P = (self.simu.K1/self.simu.M1)**0.5
        zP = self.simu.H1/(2*(self.simu.K*self.simu.M1)**0.5)
        return (w0P**2 + 2j*zP*w0P)/(w0P**2-w**2 + 2j*zP*w0P)

    def dessiner(self):
        '''Dessine le diagramme de Bode'''
        rapp = 1000
        #Reinitialisation du canvas
        for i in self.listespe:
            self.canvas.delete(i)
        l=[]
        self.w=[]
        self.f=[]
        self.f2=[]
        self.f3=[]
        if self.ideal:
            self.simu.K1 = self.simu.M1*9.81/(self.simu.L)
        #Génération des points du diagramme de Bode
        for i in range(1,int(rapp*1.1)):
            w0 = i/rapp*int(self.rapp.get())
            self.w.append(w0)
            self.f.append(20*log10(abs(self.fctDeTrsfrt(w0))))
            self.f2.append(20*log10(abs(self.fctSanspendule(w0))))
            self.f3.append(20*log10(abs(self.fctPendule(w0))))
        F = max(max(self.f), abs(min(self.f)))
        #Affichage des points du diagramme de Bode
        m=(0,0)
        m2=(0,0)
        for i in range(len(self.f)-1):
            l.append(self.canvas.create_line(20 + 560*(self.w[i]-0)/(self.w[-1]), 350 - 200*(self.f[i])/F, 20 + 560*(self.w[i+1]-0)/(self.w[-1]), 350 - 200*(self.f[i+1])/F , fill='blue'))
            if self.f[i]>m2[0]:
                m2=(self.f[i],i)
            l.append(self.canvas.create_line(20 + 560*(self.w[i]-0)/(self.w[-1]), 350 - 200*(self.f2[i])/F, 20 + 560*(self.w[i+1]-0)/(self.w[-1]), 350 - 200*(self.f2[i+1])/F , fill='purple'))
            if self.f2[i]>m[0]:
                m=(self.f2[i],i)
            l.append(self.canvas.create_line(20 + 560*(self.w[i]-0)/(self.w[-1]), 350 - 200*(self.f3[i])/F, 20 + 560*(self.w[i+1]-0)/(self.w[-1]), 350 - 200*(self.f3[i+1])/F , fill='red'))
        #Affichage des échelles w
        for i in range(4):
            l.append(self.canvas.create_text(20 +560* 10**(i-1)/self.w[-1], 360, text='w={0}'.format(10**(i-1)), fill='grey'))
            l.append(self.canvas.create_line(20 + 560*(10**(i-1))/(self.w[-1]), 0, 20 + 560*(10**(i-1)-0)/(self.w[-1]), 600, fill='grey'))
        #Affichage des échelles amplitude
        Amax = 350 - 200*(self.f2[m[1]])/F
        Amin = 350 - 200*(self.f[m2[1]])/F
        l.append(self.canvas.create_text(560, Amax-10, text='H={0}'.format(round(m[0],2)), fill='purple'))
        l.append(self.canvas.create_line(0, Amax, 600, Amax , fill='grey'))
        l.append(self.canvas.create_text(560, Amin-10, text='H={0}'.format(round(m2[0],2)), fill='blue'))
        l.append(self.canvas.create_line(0, Amin, 600, Amin , fill='grey'))
        l.append(self.canvas.create_line(0,350,600,350, width=2))
        self.listespe = l

    def modific(self):
        '''Ouvre un deuxième mode du programme Python'''
        scale = []
        try:
            self.fenetre.fenetre.destroy()
        except:
            pass
        self.fenetre._break=True
        self.bouton2.destroy()
        self.ideal = True
        # Génération de l'interface
        self.pcursor = PanedWindow(self.panedBode, orient = VERTICAL)
        self.panedCurseur = PanedWindow(self.pcursor, height=350)
        self.panedCurseur2 = PanedWindow(self.pcursor)
        s = Scale(self.panedCurseur, variable = StringVar(self.fenetreBode,self.simu.M1*1000),to=1.00,from_=1000, label='MP (en g)')
        self.panedCurseur.add(s)
        scale.append(s)
        s = Scale(self.panedCurseur, variable = StringVar(self.fenetreBode,self.simu.L*1000),to=1.00,from_=1000, label='L (en mm)')
        self.panedCurseur.add(s)
        scale.append(s)
        s = Scale(self.panedCurseur, variable = StringVar(self.fenetreBode,self.simu.H1*100),to=1.00,from_=2000, label='hP (en dg/s)')
        self.panedCurseur.add(s)
        scale.append(s)
        s = Scale(self.panedCurseur2, variable = StringVar(self.fenetreBode,self.simu.K),to=1.00,from_=1000, label='kT (en N/m)')
        self.panedCurseur2.add(s)
        scale.append(s)
        s = Scale(self.panedCurseur2, variable = StringVar(self.fenetreBode,self.simu.M*1000),to=1.00,from_=4000, label='MT (en g)')
        self.panedCurseur2.add(s)
        scale.append(s)
        s = Scale(self.panedCurseur2, variable = StringVar(self.fenetreBode,self.simu.H*1000),to=1.00,from_=2000, label='HT (en g/s)')
        self.panedCurseur2.add(s)
        scale.append(s)
        # s = Scale(self.panedCurseur1, variable = StringVar(self.fenetreBode,self.simu.K2),to=1.00,from_=1000, label='kP2 (en N/m)')
        # self.panedCurseur1.add(s)
        # scale.append(s)
        self.fenetreBode.attributes('-fullscreen',True)
        self.fenetreBode.bind('<Escape>', lambda e:self.fenetreBode.destroy())
        self.scale = scale
        self.pcursor.add(self.panedCurseur)
        self.pcursor.add(self.panedCurseur2)
        self.panedBode.add(self.pcursor)

    def after(self):
        ''' Programme qui tourne en boucle pour s'occuper de l'interface diagramme de Bode'''
        try:
            w = float(self.fenetre.W.get())
        except:
            w=0
        for i in self.liste:
            self.canvas.delete(i)
        for i in self.listeP:
            self.canvasParametres.delete(i)
        self.liste=[]
        self.listeP=[]
        if not self.ideal:
            self.liste.append(self.canvas.create_line(20 + 560*(w)/(self.w[-1]), 0, 20 + 560*(w-0)/(self.w[-1]), 600, fill='red'))
        n=10
        if self.ideal:
            self.simu.K1 = self.simu.M1*9.81/(self.simu.L)
        # Affichage des constantes dans le menu
        for i in [('ω0T',(self.simu.K/self.simu.M)**0.5, '/s'),('ω0P',((self.simu.K1+self.simu.KP)/self.simu.M1)**0.5, '/s'),('ξ',(self.simu.K1+self.simu.KP)/self.simu.K, ''),('φ',self.simu.H1/self.simu.H, ''),('ζT',self.simu.H1/(2*(self.simu.K*self.simu.M)**0.5),''),('ζP',self.simu.H1/(2*((self.simu.K1+self.simu.KP)*self.simu.M1)**0.5),''),('kT',self.simu.K,'N/m'),('kP',self.simu.K1,'N/m'),('kP2',self.simu.KP,'N/m'),('mT',self.simu.M,'kg'),('mP',self.simu.M1,'kg'),('hT',self.simu.H,'kg/s'),('hP',self.simu.H1,'kg/s'),("L",self.simu.L,'m')]:
            self.listeP.append(self.canvasParametres.create_text(75,n, text = '{0}={1} {2}'.format(i[0],round(i[1],3),i[2])))
            n+=20
        #Récupération des constantes depuis le menu défilé
        if self.ideal:
            self.simu.K1 = self.simu.M1*9.81/(self.simu.L)
            self.simu.M1 = self.scale[0].get()/1000
            self.simu.L = self.scale[1].get()/1000
            self.simu.H1 = self.scale[2].get()/100
            self.simu.K = self.scale[3].get()
            self.simu.M = self.scale[4].get()/1000
            self.simu.H = self.scale[5].get()/1000
            self.dessiner()
        self.fenetreBode.after(1,self.after)


class Interface():
    def __init__(self,th):
        ''' Interface principale'''
        self.fenetre=Tk()
        self._break=False
        self.simu=Simu(self,th)
        self.OK = False
        self.pause = True
        self.packer = False
        self.fenetre.title('Théorie et simulation de la tour : Théorie {0}'.format(th))
        # Génération de l'interface
        self.screen = PanedWindow()
        self.screen.pack(fill=BOTH, expand=1)
        self.canvas = Canvas(self.screen,height=500,width=500,bg='#A1D6E9')
        self.screen.add(self.canvas)

        self.panedDetail = PanedWindow(self.screen, orient=VERTICAL)
        self.screen.add(self.panedDetail)
        lab = Label(self.panedDetail, text="Taille pendule (en m)")
        self.panedDetail.add(lab)
        self.tPendule = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.tPendule)
        self.panedDetail.add(e)
        self.tPendule.set("0.25")
        if self.simu.theorie == 0:
            lab = Label(self.panedDetail, text="X (en m)")
            self.panedDetail.add(lab)
        if self.simu.theorie == 1:
            lab = Label(self.panedDetail, text="Xt (en m)")
            self.panedDetail.add(lab)
        self.tx = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.tx)
        self.panedDetail.add(e)
        self.tx.set("0")
        if self.simu.theorie == 0:
            lab = Label(self.panedDetail, text="U (en rad)")
            self.panedDetail.add(lab)
        if self.simu.theorie == 1:
            lab = Label(self.panedDetail, text="Xp (en m)")
            self.panedDetail.add(lab)
        self.tu = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.tu)
        self.panedDetail.add(e)
        self.tu.set("0")

        lab = Label(self.panedDetail, text="mT : Masse tour (en kg)")
        self.panedDetail.add(lab)

        self.m = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.m)
        self.panedDetail.add(e)
        self.m.set("1.5")

        lab = Label(self.panedDetail, text="mP : Masse du pendule (en kg)")
        self.panedDetail.add(lab)

        self.m1 = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.m1)
        self.panedDetail.add(e)
        self.m1.set("0.3")

        lab = Label(self.panedDetail, text="hT : Frottement tour (en kg/s)")
        self.panedDetail.add(lab)

        self.h = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.h)
        self.panedDetail.add(e)
        self.h.set("0.0857")

        lab = Label(self.panedDetail, text="hP : Frottements pendule (en kg/s)")
        self.panedDetail.add(lab)

        self.h1 = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.h1)
        self.panedDetail.add(e)
        self.h1.set("10.8")

        lab = Label(self.panedDetail, text="kT : Rappel élastique tour (en N/m)")
        self.panedDetail.add(lab)

        self.k = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.k)
        self.panedDetail.add(e)
        self.k.set("99")

        self.k1 = StringVar(self.panedDetail)
        if self.simu.theorie == 0:
            lab = Label(self.panedDetail, text="kP : Rappel élastique pendule (en N/m)")
            self.panedDetail.add(lab)
            e = Entry(self.panedDetail, textvariable = self.k1)
            self.panedDetail.add(e)
        self.k1.set("5")

        self.kP = StringVar(self.panedDetail)
        lab = Label(self.panedDetail, text="kP2 : Rappel élastique supplémentaire (en N/m)")
        self.panedDetail.add(lab)
        e = Entry(self.panedDetail, textvariable = self.kP)
        self.panedDetail.add(e)
        self.kP.set("0")

        self.A = StringVar(self.panedDetail)
        if self.simu.theorie==0:
            lab = Label(self.panedDetail, text="Accélération du sol (en m/s²)")
            self.panedDetail.add(lab)
            self.A.set("2")
        if self.simu.theorie==1:
            lab = Label(self.panedDetail, text="Amplitude du sol (en m)")
            self.panedDetail.add(lab)
            self.A.set("0.001")

        e = Entry(self.panedDetail, textvariable = self.A)
        self.panedDetail.add(e)

        lab = Label(self.panedDetail, text="Pulsation forcée (en /s)")
        self.panedDetail.add(lab)

        self.W = StringVar(self.panedDetail)
        e = Entry(self.panedDetail, textvariable = self.W)
        self.panedDetail.add(e)
        if self.simu.theorie == 0:
            self.W.set("0")
        if self.simu.theorie == 1:
            self.W.set("10")

        self.bouton = Button(self.panedDetail,text='Valider', fg = 'green', command=self.valid)
        self.panedDetail.add(self.bouton)

        self.panedPause = PanedWindow(self.panedDetail, orient=HORIZONTAL)
        self.panedDetail.add(self.panedPause)

        self.boutonmoins = Button(self.panedPause,text=' - Temps ', fg = 'blue', command=self.mettremoins)
        self.panedPause.add(self.boutonmoins)

        self.boutonpause = Button(self.panedPause,text='Pause/Play', fg = 'red', command=self.mettrepause)
        self.panedPause.add(self.boutonpause)

        self.boutonplus = Button(self.panedPause,text=' + Temps ', fg = 'blue', command=self.mettreplus)
        self.panedPause.add(self.boutonplus)
        # Génération de l'interface qui affiche les courbes
        self.screenCourbe = PanedWindow(self.screen, orient=VERTICAL)
        self.canvascourbe = Canvas(self.screenCourbe, width=500,height=300, bg = 'white')
        self.screenCourbe.add(self.canvascourbe)
        self.canvascourbependule = Canvas(self.screenCourbe, width=500, height = 200, bg='white')
        self.screenCourbe.add(self.canvascourbependule)
        self.labelMin = Label(self.screenCourbe, text = 'Minimum')
        self.screenCourbe.add(self.labelMin)
        self.tempsmin = StringVar(self.screen,'auto')
        self.entreMin = Entry(self.screenCourbe, text=self.tempsmin)
        self.screenCourbe.add(self.entreMin)
        self.labelMax = Label(self.screenCourbe, text = 'Maximum')
        self.screenCourbe.add(self.labelMax)
        self.tempsmax = StringVar(self.screen,'auto')
        self.entreMax = Entry(self.screenCourbe, text=self.tempsmax)
        self.screenCourbe.add(self.entreMax)
        self.labelDefaut = Label(self.screenCourbe, text = 'Défaut (meilleure qualité de courbe à 1)')
        self.screenCourbe.add(self.labelDefaut)
        self.defaut = StringVar(self.screen,'10')
        self.defaute = Entry(self.screenCourbe, text=self.defaut)
        self.screenCourbe.add(self.defaute)
        self.boutonCech = Button(self.screenCourbe, text = 'Changer échelle', fg='orange', command = self.modifEchelle)
        self.screenCourbe.add(self.boutonCech)
        self.boutonBode = Button(self.screenCourbe, text = 'Lancer Interface Bode', command = self.launchbode, fg='purple')
        self.screenCourbe.add(self.boutonBode)
        self.screen.add(self.screenCourbe)
        # Génération des variables initiales nécessaires au programme
        self.t0 = 0 #temps de la simulation
        self.t1 = 0
        self.AV = 0.02 #avancement temporel
        self.N=10
        self.t2 = None #compteur temporel de la tour 1
        self.t3 = None #compteur temporel de la tour 2
        self.list_coord = {0:[],1:[],2:[]}
        self.list_coordpendule = []
        self.list_obj = []
        self.modifEchelle()
        self.valid()
        self.k1.set(str(self.simu.M1*9.81/self.simu.L))
        self.after()
        self.fenetre.attributes('-fullscreen',True)
        self.fenetre.bind('<Escape>', lambda e:self.fenetre.destroy())
        self.fenetre.mainloop()

    def modifEchelle(self):
        '''Commande pour modifier l'échelle du graphe'''
        self.N = int(self.defaut.get())+1
        self.tmin=0
        if self.tempsmin.get() != 'auto':
            self.tmin=int(self.tempsmin.get())
        if self.tempsmax.get()=='max' or self.tempsmax.get()=='auto':
            self.tmax = 10000
            if self.tempsmin.get() == 'auto' and self.list_coord[2] != []:
                self.tmin = int(self.list_coord[2][-1]-10)
        else:
            self.tmax = int(self.tempsmax.get())
            if self.tempsmin.get() == 'auto':
                self.tmin = self.tmax-10

        if self.tmax <= self.tmin:
            self.tempsmax.set(str(self.tmin+10))
            self.tmax = self.tmin + 10
        if self.list_coord[2] != []:
            self.packer = False

    def launchbode(self):
        '''Commande pour lancer l'interface du diagramme de Bode'''
        Bode(self)

    def valid(self):
        '''lorsque l'utilisateur clique sur play/pause'''
        try:
            # Réinitialisation de la simulation
            self.simu.Eul.X=[0,0,0]
            self.simu.Eul.U=[0,0,0]
            self.simu.EulS.X=[0,0,0]
            self.simu.Eul.X[0] = float(self.tx.get())
            self.simu.EulS.X[0] = float(self.tx.get())
            self.simu.Eul.U[0] = float(self.tu.get())
            self.simu.L = float(self.tPendule.get())
            self.simu.M = float(self.m.get())
            self.simu.M1 = float(self.m1.get())
            self.simu.H = float(self.h.get())
            self.simu.H1 = float(self.h1.get())
            self.simu.K = float(self.k.get())
            if self.simu.theorie == 0:
                self.simu.K1 = float(self.k1.get())
            if self.simu.theorie == 1:
                self.simu.K1 = self.simu.M1*9.81*(0.3**2+self.simu.Eul.U[0])/(0.3**2)
            self.simu.KP = float(self.kP.get())
            a,b = float(self.A.get()),float(self.W.get())
            self.simu.Fex = lambda t:a* sin(b*t)
            self.t1 = 0
            self.t2 = None
            self.t3 = None
        except:
            # la raison pour laquelle il pourrait y avoir un problème est si l'utilisateur ne rentre pas un nombre dans les entrées exemple ' '
            pass
        # Réinitialisation des coordonnées
        for i in [0,1,2]:
            self.list_coord[i]=[]
        self.list_coordpendule=[]
        self.pause = False
        self.t0=0

    def mettrepause(self):
        '''Commande du bouton pause'''
        self.pause = not self.pause

    def mettremoins(self):
        '''Réduire le temps de défilement de l'interface'''
        if self.AV < 0.001:
            self.AV = 0.001
        self.AV *= 0.5

    def mettreplus(self):
        '''Augmente le temps de défilement de l'interface'''
        self.AV *= 2

    def after(self):
        '''Boucle qui fait tourner le programme en permanence'''
        AVANCEMENTtps = self.AV
        # quand on mets pause il faut réinitialiser la courbe
        if self.tempsmin.get() == 'auto' and self.list_coord[2] != []:
            self.tmin = int(self.list_coord[2][-1]-10)
        if self.pause and not self.packer and self.list_coord[2]!=[]:
            self.packer = True
            # on commence par supprimer l'ancienne courbe :
            for i in self.list_obj:
                self.canvascourbe.delete(i)
                self.canvascourbependule.delete(i)
            self.list_obj=[]

            # on trace la nouvelle courbe
            h1=self.canvascourbe.create_line(480,280,20,280)
            h2=self.canvascourbe.create_line(20,20,20,280)
            T=max(self.list_coord[2]) #maximum d'amplitude tour
            Hmax = max([max(self.list_coordpendule), abs(min(self.list_coordpendule))]) #maximum amplitude pendule
            a=min(self.tmax,self.list_coord[2][-1])
            if a<self.tmin:
                a=self.tmin
            h3=self.canvascourbe.create_text(480,290, text=str(round(a,2))+' s')
            h4 = self.canvascourbe.create_text(20, 290, text = str(self.tmin)+' s')
            h5 = self.canvascourbe.create_text(450, 20, text = 'tour avec pendule', fill = 'blue')
            h6 = self.canvascourbe.create_text(450, 40, text = 'tour sans pendule', fill = 'purple')
            # on prend le maximum d'amplitude
            X=max(max(max(self.list_coord[0]),abs(min(self.list_coord[0]))), max(max(self.list_coord[1]), abs(min(self.list_coord[1]))))
            h7 = self.canvascourbe.create_text(40, 10, text = '{0} mm'.format(round(X*10**(3),3)))
            h8 = self.canvascourbependule.create_line(10,190,490,190)
            h9 = self.canvascourbependule.create_line(10,10,10,190)
            h10 = self.canvascourbependule.create_text(60, 10, text = 'max à '+str(round(Hmax*10**(3),3))+' mm')
            h11 = self.canvascourbependule.create_text(450, 20, text = 'pendule', fill = 'red')
            h13=self.canvascourbependule.create_text(480,190, text=str(round(a,2))+' s')
            h14 = self.canvascourbependule.create_text(20, 190, text = str(self.tmin)+' s')
            if X!= 0:
                N = self.N
                tmin,tmax = self.tmin,self.tmax
                # par soucis de complexité et d'affichage, on choisit le rapport de point pris (1/N)
                N1,N2 = recherche(self.list_coord[2],tmin,tmax)
                taille = N2-N1
                for j in range(len(self.list_coord[0])-N):
                    if j%N == 0:
                        a,b=self.list_coord[0][j],self.list_coord[0][j+N]
                        c,d=self.list_coord[1][j],self.list_coord[1][j+N]
                        e,f = self.list_coordpendule[j], self.list_coordpendule[j+N]
                        t = self.list_coord[2][j]
                        if tmin <= t and t <= tmax:
                            hj = self.canvascourbe.create_line(20+(j-N1)*460/taille, 150-130*a/X, 20+(j-N1+N)*460/taille, 150-130*b/X, fill='blue')
                            hj2 = self.canvascourbe.create_line(20+(j-N1)*460/taille, 150-130*c/X, 20+(j-N1+N)*460/taille, 150-130*d/X, fill='purple')
                            hj3 = self.canvascourbependule.create_line(20+(j-N1)*460/taille, 100-80*e/Hmax, 20+(j-N1+N)*460/taille, 100-80*f/Hmax, fill='red')
                            self.list_obj.append(hj)
                            self.list_obj.append(hj2)
                            self.list_obj.append(hj3)
            for i in [h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h13,h14]:
                 self.list_obj.append(i)

        if not self.pause:
            # si on n'est pas en pause, on n'affiche pas la courbe et on modifie la valeurs des points
            self.packer = False
            self.t1 += AVANCEMENTtps
            # On récupère les nouvelles coordonnées à l'aide de la méthode d'Euler
            T=self.t0
            H=0.02/(2**6)
            while T-self.t0 < AVANCEMENTtps:
                T+=H
                self.simu.Eul.PointSuivant(T,H)
                self.simu.EulS.PointSuivant(T,H)
                self.list_coord[0].append(self.simu.Eul.X[0])
                self.list_coord[1].append(self.simu.EulS.X[0])
                self.list_coord[2].append(T)
                self.list_coordpendule.append(self.simu.Eul.U[0])
            self.t0 = T

            # On modifie le mouvement de la tour
            self.simu.tourS.pack()
            self.simu.tour.pack()

            try:
                self.canvas.delete(self.labelTemps)
                self.canvas.delete(self.labelTemps2)
                self.canvas.delete(self.labelAmplitude)
                self.canvas.delete(self.labelAmplitudeS)
            except:
                pass

            # on modifie les compteurs temporels
            if self.t2 == None:
                self.labelTemps = self.canvas.create_text(125,50, text=str(round(self.t1,4))+' s')
                self.labelAmplitude = self.canvas.create_text(125,70, text='X = ' +str(round(self.simu.Eul.X[0],4))+' m')
            else :
                self.labelTemps = self.canvas.create_text(125,50, text=str(round(self.t2,4))+' s')

            if self.t3 == None:
                self.labelTemps2 = self.canvas.create_text(375,50, text=str(round(self.t1,4))+' s')
                self.labelAmplitudeS = self.canvas.create_text(375, 70, text='X = '+str(round(self.simu.EulS.X[0],4))+' m')
            else :
                self.labelTemps2 = self.canvas.create_text(375,50, text=str(round(self.t3,4))+' s')

            if self.t2 != None and self.t3 != None:
                self.pause = True


            self.simu.tour.x = 400*x(self.simu.Eul.X)
            self.simu.tour.u = 400*x(self.simu.Eul.U)
            self.simu.tourS.x = 400*x(self.simu.EulS.X)

            # moment où le compteur s'arrête
            PRECISION = 5
            if abs(self.simu.tour.x) < 10**(-PRECISION) and abs(self.simu.tour.u) < 10**(-PRECISION) and abs(self.simu.Eul.X[1]) < 10**(-PRECISION) :
                self.simu.tour.x = 0
                self.simu.tour.u = 0
                self.simu.Eul.X,self.simu.Eul.U = [0,0,0], [0,0,0]
                if self.t2 == None:
                    self.t2 = self.t1

            if abs(self.simu.tourS.x) < 10**(-PRECISION) and abs(self.simu.EulS.X[1]) < 10**(-PRECISION) :
                self.simu.tourS.x = 0
                self.simu.tourS.u = 0
                self.simu.EulS.X,self.simu.EulS.U = [0,0,0], [0,0,0]
                if self.t3 == None:
                    self.t3 = self.t1
        if not self._break:
            self.fenetre.after(1,self.after)

def recherche(liste,i1,i2):
    k1,k2=0,0
    for i in range(len(liste)) :
        if liste[i]<i1:
            k1=i
        if liste[i]<i2:
            k2=i
        else:
            return k1,k2
    return k1,k2+0.01 #pour pas diviser par zéro

class StartInterface():
    def __init__(self):
        '''Interface de départ'''
        # Cette interface propose deux théories différentes
        self.fenetre = Tk()
        self.fenetre.title('Quelle simulation choisir ?')
        self.paned=PanedWindow(self.fenetre, orient = VERTICAL)
        self.bouton0 = Button(self.paned, text = "Théorie 0", command = self.theorie0)
        self.paned.add(self.bouton0)
        self.lab = Label(self.paned, text = 'Modèle du pendule simple : oscillation libre proche de la réalité, prédiction des oscillations forcées pas convenable', bg='cyan')
        self.paned.add(self.lab)
        self.bouton1 = Button(self.paned, text = "Théorie 1", command = self.theorie1)
        self.paned.add(self.bouton1)
        self.lab = Label(self.paned, text = 'Modèle du TMD classique, diagramme de Bode idéal', bg='lightgreen')
        self.paned.add(self.lab)
        self.paned.pack()
        self.fenetre.mainloop()

    def theorie0(self):
        '''L'utilisateur a choisi la théorie 0'''
        self.fenetre.destroy()
        Interface(0)

    def theorie1(self):
        '''L'utilisateur a choisi la théorie 1'''
        self.fenetre.destroy()
        Interface(1)

StartInterface()