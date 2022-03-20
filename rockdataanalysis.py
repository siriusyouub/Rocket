import math 
import numpy as np
import random
M_L = 1 # kg payload mass
list_h_max = [float(10000*0.3048), float(20000*0.3048), float(30000*0.3048)] # maximum altitude, unit m
list_a_max = [5,10, 20] # normalized maximum acceleration
list_R_max = [6,11,21] # R = a+1
list_SM = [1,2,3] # static margin = (X_cp - X_cg)/D, x_cp--distance of center of pressure to the head, x_cg--distance of center of gravity to the head
rho_s = 2700 # kg/m^3 shell density (aluminum)
rho_p = 1772 # propellant density, 66 - 78% AP, 18% organic polymer, 4 - 20% Al
sigma_s = 60000000 # Pa shell working stress
N = 3 # number of fins
g = 9.81
pi = float(3.14)
#R =  Mo/Mb # mass ratio R = 1+M_p/(M_s+M_L)
a = float(math.sqrt(1.4*287*298))
#W_eq = I_sp*g # v_eq = I_sp*t_b
#Cp = (P-P_inf)/(0.5*rho_inf*(v_inf**2))
#F_p = 0.5*rho_inf*Cp*Ap*(v_inf**2) # Cp is press coeff 
#W_eq = math.sqrt((h_max*g)/((ln(R_max-1)/2)*(ln(R_max)-2)+(R_max-1)/(R_max))) # 这里的 R_max 和 h_max 是列表，需要调用一下,算完后需要给它新建一个列表

Pa = 101325 # Pa
#fin = random.uniform(start = 10, end= 15) #单位m, for D

i = 0
a = 0
b = 0
c =0
while True: 

    for a in range(3):
        h_max = list_h_max[a]
        for b in range(3):
            a_max = list_a_max[b]
            R_max = list_R_max[b]
            for c in range(3):
                SM = list_SM[c]
                Weq = math.sqrt((h_max*g)/(0.5*math.log(R_max)*(math.log(R_max)-2)+a_max/R_max))
                Me = Weq/a_max
                P0 = Pa*(1+0.2*(Me**2))**(1.4/0.4)
                #P0 = 1.26*Pa
                d = 0
                while True:
                    
                    list_D = np.linspace((0.01),(0.1),int(450))
                    #for D in (np.linspace((0.01),(0.1),int(450))):

                    while True: 
                        e = 0
                        list_L = np.linspace((0.6),(0.95),int(2501))
                        D = list_D[d]

                        while True:
                            list_L_y = list()
                            list_Ms = list()
                            list_Mp = list()
                            L = list_L[e]                            
                            #for L in np.linspace((0.6),0.95,int(2500)):
                            delta = D * P0/(2*sigma_s)

                            M_s = rho_s*((D+L)*(pi/4)*delta*(D**2)+pi*(D/2)*(math.sqrt(5/4)*D)*delta+3*delta*(D**2)/2) #Asumme the rocket is thin surfaace M_s = F(D,L,delta,rho_s) M_s = M_o-M_p-M_L

                            M_p = float((R_max-1)*(M_s+M_L))
                            L_p = M_p/(pi*(D**2)*rho_p/4)
                            y = L + D-L_p
                            if y > 0 & e<=2500:
                                list_L_y.append(L)
                                list_Mp.append(M_p)
                                list_Ms.append(M_s)
                                e = e+1
                                continue
                            elif y<=0 & e<=1500:
                                e = e+1
                                continue
                            elif e>2500:
                                break
                        f = 0
                        #for f in range(list_L_y):
                        while True:
                            list_XcgXxpDSM = list()
                            L_y = list_L_y[f]
                            Mp = list_Mp[f]
                            Ms = list_Ms[f]
                            list_MP = list()
                            list_MS = list()

                            A2 = float(D*pi*L_y) #surface area for body
                            A3 = float(D*pi*D) #surface area for bottom
                            #A3 = float()

                            Ap1 = float(0.25*(D**2)) # project area of nose
                            Ap2 = float(D*L_y) #project area of body
                            Ap3 = float(2*(D**2)) #project area of bottom

                            C_Nf = 18.5 #used the table in page 16 of calculate center of pressure
                            S = D # length of fin
                            R_rocket = D/2 # radius of rocket body 
                            K_fb = 1+(R_rocket/S)/(1+R_rocket/S)

                            Cp1 = 2 # nose
                            Cp2 = 0 # body
                            Cp3 = float(12/(1+math.sqrt(6)))#use this function to get float(4*n*((s/d)**2)/(1+math.sqrt(1+(2*l/(a+b))**2))) 
                            #bottom, used the table in page 16 of calculate center of pressure
                            # it is a constant value


                            M1 = float((math.sqrt(5)/2)*D*pi*(D/2)*rho_s*delta) # mass of nose only shell
                            M2 = float(rho_s*L*pi*D*delta) #mass of body shell
                            M3 = float(rho_s*pi*delta*(D**2)) # mass of bottom, include 3 fins
                            M4 = float(3*delta*rho_s*(D**2)/2)#fin
                            M5 = 1#payload
                            M6 = float(L_p*rho_p*(1/4)*pi*D**2)#propollent

                            X_cp1 = float((2/3)*D) # center of pressure for nose
                            X_cp2 = float(0.5*L+0.25*M2*(D**2)/(A2*L))
                            X_cp3 = float(L+D+0.5*D) #0.5 get from chart 6, use the project page 26
                            #X_cp2 = float(0.5*L+0.25*D*(M_s+(M_p*(L_p-D)/L_p))/(pi*L**2))# center of pressure for body
                            #X_cp3 = float(0.5*D+0.25*(M_s+D*M_p/L_p)/(D*pi))

                            X_cg1 = float((2/3)*D) #center of gravity for nose
                            X_cg2 = float(D+L/2) #center of gravity for body
                            X_cg3 = float(D+L+D/2) #center of gravity for bottom
                            X_cg4 = float(L_p-D/3)#fin
                            X_cg5 = float((2/3)*D)#payload
                            X_cg6 = float(2*D+L-L_p/2)#propollent

                            #X_cg = [] # central gravity for seperate part, need list
                            X_cps = float((X_cp1*Cp1+X_cp3*Cp3)/(Cp1+Cp3))#F(D,L) 
                            #X_cps = float((X_cp1*Cp1*Ap1+X_cp3*Cp3*Ap3)/(Cp1*Ap1+Cp3*Ap3))#F(D,L)
                            X_cgs = float((X_cg1*M1+X_cg2*M2+X_cg3*M3+X_cg4*M4+X_cg5*M5+X_cg6*M6)/(M_p))#F(D,L,L_p,rho_p,rho_s,delta_s), total X_cg

                            R= 1+M_p/(M_s+M_L)
                            tb = (R_max-1)*Weq/(g*R_max) #需要建一个空白的列表，调用R_max
                            #print(X_cps-X_cgs-D*SM)
                            #Xcg_X_xp_DSM = abs(X_cps-X_cgs-D*SM)
                            #print(Xcg_X_xp_DSM)
                        if (X_cps-X_cgs-D*SM) < 0.0001 & f<= range(list_L_y): #use the function mimimum
                            f = f+1
                            list_MS.append(Ms)
                            list_MP.append(Mp)

                            continue
                        elif (X_cps-X_cgs-D*SM) >=0.0001 & f<= range(list_L_y): 
                            f = f+1
                            continue
                        elif f> range(list_L_y):
                            break
                    MS = list_MS
                    MP = list_MP
                    g = 0
                    while g<=range(list_MS):

                        list_lamda = float(M_L/(MS+MP))
                        if list_lamda[g] != max(list_lamda) & g<=(range(list_MS)):
                            g = g+1
                            continue
                        elif list_lamda[g] != max(list_lamda) & g>(range(list_MS)):
                            
                            break 

                        elif list_lamda[g] == max(list_lamda):
                            lamda = list_lamda[g]
                            print('h_max='+str(list_h_max[a])+'m', end=" ")
                            print('a_max='+str(list_a_max[b])+'m', end=" ")
                            print('SM='+str(list_SM[c])+'m', end=" ")
                            print("""a_max={}\n SM = {}\n R = {}\n Weq = {}\n tb = {}\n P0/Pa={}\n delta/D={}\n D={}\n L/D={}\n X_cgs={}\n X_cps={}\n M_p={}\n M_s={}\n M_o={}\n lamda={}\n Epsilon={}\n""".format(a_max, SM, R, Weq, tb, P0/Pa, delta/D, D, L/D, X_cgs,X_cps,M_p,M_s,M_o,namda,E))
                        #print()
                            break

                        
                

    
    #print('R = '+str(R)+ 'Weq = '+str(Weq)+ 'tb = '+ str(tb),'P0 = '+str(P0)+'delta = '+ str(delta/D),
    #'D='+str(D),'L = '+str(L),'X_cgs='+str(X_cgs)+'X_cps='+str(X_cps)+'M_p'+str(M_p)+'M_o'+str(M_o)+'lamda'+str(namda)+'E='+str(E))

    #print("""R ={}, W_eq={}, t_b={}, P0={}, delta/D={}, D={}
        # , L={}, X_cg={}, X_cp={},
       # M_p={},M_o={},namda={}""".format(str(R, W_eq, tb,P0,delta/D,D,L,X_cgs,X_cps,M_p,M_o,namda,E)))