'''
Damping coefficients components acc. to IKEDA method with further improvements

Kawahara, Y., Maekawa, K., Ikeda, Y. (2011). A Simple Prediction Formula of Roll Damping of Conventional Cargo Ships on the Basis of Ikeda’s Method and Its Limitation. In: Almeida Santos Neves, M., Belenky, V., de Kat, J., Spyrou, K., Umeda, N. (eds) Contemporary Ideas on Ship Stability and Capsizing in Waves. Fluid Mechanics and Its Applications, vol 97. Springer, Dordrecht. https://doi.org/10.1007/978-94-007-1482-3_26

'''
import numpy as np

def bilke(cb,cm,ogd,brth,bbk,lbk,lpp,t,phi,tw):
    G    = 9.80665
    bd   = brth/t
    lbkl = lbk/lpp
    bbkb = bbk/brth
    om   = 2*np.pi/tw
    omd  = om*np.sqrt(brth/2/G)

    fbk1  = (-0.3651*cb+0.3907)*(bd-2.83)**2-2.21*cb+2.632
    fbk2  = 0.00255*phi**2+0.122*phi+0.4794
    fbk3  = (-0.8913*bbkb**2-0.0733*bbkb)*lbkl**2+(5.2857*bbkb**2-0.01185*bbkb+0.00189)*lbkl
    abk   = fbk1*fbk2*fbk3
    bbk1  = (5.0*bbkb+0.3*bd-0.2*lbkl+0.00125*phi**2-0.0425*phi-1.86)*ogd
    bbk2  = -15.0*bbkb+1.2*cb-0.1*bd-0.0657*ogd**2+0.0586*ogd+1.6164
    bbk3  = 2.5*ogd+15.75
    return abk*np.exp(bbk1+bbk2*cm**bbk3)*omd

def eddy(cb,cm,ogd,brth,t,phi,tw):
    G   = 9.80665
    phi = np.radians(phi)
    om  = 2*np.pi/tw
    omd = om*np.sqrt(brth/2/G)
    bd  = brth/t

    fe1  = (-0.0182*cb+0.0155)*(bd-1.8)**3
    fe2  = -79.414*cb**4+215.695*cb**3-215.883*cb**2+93.894*cb-14.848
    ae   = fe1 + fe2
    be1  = (3.98*cb-5.1525)*(-0.2*bd+1.6)*ogd*((0.9717*cb**2-1.55*cb+0.723)*ogd+0.04567*cb+0.9408)
    be2  = (0.25*ogd+0.95)*ogd-219.2*cb**3+443.7*cb**2-283.3*cb+59.6
    be3  = -15*cb*bd+46.5*cb+11.2*bd-28.6
    cr   = ae*np.exp(be1+be2*cm**be3)
    return 4*omd*phi/(3*np.pi*cb*bd**3)*cr

def fric(lpp,brth,t,cb,ogd,phi,tw):
    G   = 9.80665
    phi = np.radians(phi)
    rho = 102.59
    nu  = 1.14e-6
    om  = 2*np.pi/tw
    bd  = brth/t
    rf  = t*((0.887+0.145*cb)*(1.7+cb*bd)-2.0*ogd)/np.pi
    sf  = lpp*(1.75*t+cb*brth)
    cf  = 1.328*((3.22*rf**2*phi**2)/(tw*nu))**(-.5)
    bf  = 4/3*np.pi*rho*sf*rf**3*phi*om*cf
    return bf/(rho*lpp*brth**3*t*cb)*np.sqrt(brth/2/G)

def lift(lpp,brth,t,cb,cm,ogd,u):
    G   = 9.80665
    rho = 1025.9
    nu  = 1.14e-6
    bd  = brth/t
    if cm<=.92:
        kp = 0.0
    elif (cm<=.97) & (cm>.92):
        kp = 0.1
    else:
        kp = 0.3

    kn  = 2*np.pi*t/lpp+kp*(4.1*brth/lpp+0.045)
    og  = ogd*t
    lo  = 0.3*t
    lr  = 0.5*t
    bl  = 0.5*rho*u*lpp*t*kn*lo*lr*(1-1.4*og/lr+0.7*og**2/lo/lr)
    return bl/(rho*lpp*brth**3*t*cb)*np.sqrt(brth/2/G)

def wave(cb,cm,ogd,brth,t,phi,tw,u):
    G   = 9.80665
    phi = np.radians(phi)
    om  = 2*np.pi/tw
    bd  = brth/t

    X1,X2,X3,X4,X5 = bd,cb,cm,1-ogd,om*np.sqrt(brth/2/G)

    A111=-0.002222*X1**3+0.040871*X1**2-0.286866*X1+0.599424
    A112=0.010185*X1**3-0.161176*X1**2+0.904989*X1-1.641389
    A113=-0.015422*X1**3+0.220371*X1**2-1.084987*X1+1.834167
    A121=-0.0628667*X1**4+0.4989259*X1**3+0.52735*X1**2-10.7918672*X1+16.616327
    A122=0.1140667*X1**4-0.8108963*X1**3-2.2186833*X1**2+25.1269741*X1-37.7729778
    A123=-0.0589333*X1**4+0.2639704*X1**3+3.1949667*X1**2-21.8126569*X1+31.4113508
    A124=0.0107667*X1**4+0.0018704*X1**3-1.2494083*X1**2+6.9427931*X1-10.2018992
    A131=0.192207*X1**3-2.787462*X1**2+12.507855*X1-14.764856
    A132=-0.350563*X1**3+5.222348*X1**2-23.974852*X1+29.007851
    A133=0.237096*X1**3-3.535062*X1**2+16.368376*X1-20.539908
    A134=-0.067119*X1**3+0.966362*X1**2-4.407535*X1+5.894703

    A11=A111*X2**2+A112*X2+A113
    A12=A121*X2**3+A122*X2**2+A123*X2+A124
    A13=A131*X2**3+A132*X2**2+A133*X2+A134

    AA111=17.945*X1**3-166.294*X1**2+489.799*X1-493.142
    AA112=-25.507*X1**3+236.275*X1**2-698.683*X1+701.494
    AA113=9.077*X1**3-84.332*X1**2+249.983*X1-250.787
    AA121=-16.872*X1**3+156.399*X1**2-460.689*X1+463.848
    AA122=24.015*X1**3-222.507*X1**2+658.027*X1-660.665
    AA123=-8.56*X1**3+79.549*X1**2-235.827*X1+236.579

    AA11=AA111*X2**2+AA112*X2+AA113
    AA12=AA121*X2**2+AA122*X2+AA123

    AA1=(AA11*X3+AA12)*(1-X4)+1.0

    A1=(A11*X4**2+A12*X4+A13)*AA1
    A2=-1.402*X4**3+7.189*X4**2-10.993*X4+9.45

    A31=-7686.0287*X2**6+30131.5678*X2**5-49048.9664*X2**4+42480.7709*X2**3-20665.147*X2**2+5355.2035*X2-577.8827
    A32=61639.9103*X2**6-241201.0598*X2**5+392579.5937*X2**4-340629.4699*X2**3+166348.6917*X2**2-43358.7938*X2+4714.7918
    A33=-130677.4903*X2**6+507996.2604*X2**5-826728.7127*X2**4+722677.104*X2**3-358360.7392*X2**2+95501.4948*X2-10682.8619
    A34=-110034.6584*X2**6+446051.22*X2**5-724186.4643*X2**4+599411.9264*X2**3-264294.7189*X2**2+58039.7328*X2-4774.6414
    A35=709672.0656*X2**6-2803850.2395*X2**5+4553780.5017*X2**4-3888378.9905*X2**3+1839829.259*X2**2-457313.6939*X2+46600.823
    A36=-822735.9289*X2**6+3238899.7308*X2**5-5256636.5472*X2**4+4500543.147*X2**3-2143487.3508*X2**2+538548.1194*X2-55751.1528
    A37=299122.8727*X2**6-1175773.1606*X2**5+1907356.1357*X2**4-1634256.8172*X2**3+780020.9393*X2**2-196679.7143*X2+20467.0904

    AA311=(-17.102*X2**3+41.495*X2**2-33.234*X2+8.8007)*X4+36.566*X2**3-89.203*X2**2+71.8*X2-18.108

    AA31=(-0.3767*X1**3+3.39*X1**2-10.356*X1+11.588)*AA311
    AA32=-0.0727*X1**2+0.7*X1-1.2818

    XX4=X4-AA32

    AA3=AA31*(-1.05584*XX4**9+12.688*XX4**8-63.70534*XX4**7+172.84571*XX4**6-274.05701*XX4**5+257.68705*XX4**4-141.40915*XX4**3+44.13177*XX4**2-7.1654*XX4-0.0495*X1**2+0.4518*X1-0.61655)

    A3=A31*X4**6+A32*X4**5+A33*X4**4+A34*X4**3+A35*X4**2+A36*X4+A37+AA3


    bwdl = A1/X5*np.exp(-A2*(np.log(X5)-A3)**2/1.44)

    if u >0:
        zd  = om**2*t/G
        tau = u*om/G
        A1  = 1.0+zd**(-1.2)*np.exp(-2*zd)
        A2  = 0.5+zd**(-1.0)*np.exp(-2*zd)
        cbw = 0.5*(((A2+1)+(A2-1)*np.tanh(20*(tau-0.3)))+(2*A1-A2-1)*np.exp(-150*(tau-0.25)**2))
    else:
        cbw = 1.0

    return bwdl*cbw
