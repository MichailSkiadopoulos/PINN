from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np
import random

random.seed(6742 * 1)
np.random.seed(6742 * 1)

#------------------------------------------------------------------------------
###############################################################################
def CheckOverlap(x,y,z,r,elemsize,xi,yi,zi,ri):
    ""
    nc=len(xi)
    gamma=1
    B=True
    
    for i in range(0,nc):
        if np.sqrt((x-xi[i])**2+(y-yi[i])**2+(z-zi[i])**2)<r+ri[i] + gamma * elemsize:
            B=False
            break
    
    return B;

def CheckBorder(x,y,z,r,elemsize, W_connecting, H_connecting, L_connecting,W,H,L):
    ""
    dist1=np.array([y,z,H-y,L-z])
    dist2=np.array([x,W-x])
    gamma1=1
    gamma2=1
    gamma3=1
    A=True
    
    if abs(x)<W_connecting+gamma1*elemsize + r or abs(W-x)<W_connecting+gamma1*elemsize + r \
        or abs(y)<H_connecting+gamma2*elemsize + r or abs(H-y)<H_connecting+gamma2*elemsize + r \
        or abs(z)<L_connecting+gamma3*elemsize + r or abs(L-z)<L_connecting+gamma3*elemsize + r:
        A=False

    return A;
###############################################################################

#------------------------------------------------------------------------------
###############################################################################
def Faces_detection(H_sample, W_sample, W_connecting1, W_inf, L_sample):

    a = mdb.models['Model-%d' %(q)].rootAssembly                
    f = a.instances['Part-1-1'].faces
              
    Face_Left1 = f.findAt((W_sample + W_connecting1 / 2, H_sample / 2, L_sample))
    Face_Right1 = f.findAt((W_sample + W_connecting1 / 2, H_sample / 2, 0.0))
    Face_Left2 = f.findAt((W_sample + W_connecting1  + W_inf / 2, H_sample / 2, L_sample))
    Face_Right2 = f.findAt((W_sample + W_connecting1  + W_inf / 2, H_sample / 2, 0.0))
    
    q1 = Face_Left1.index
    q2 = Face_Right1.index
    q3 = Face_Left2.index
    q4 = Face_Right2.index
    
    F_Left1 = f[q1: q1 + 1]
    F_Right1= f[q2: q2 + 1]
    F_Left2 = f[q3: q3 + 1]
    F_Right2 = f[q4: q4 + 1] 
    
    a.Set(faces=F_Left1 , name='Face_Left1')
    a.Set(faces=F_Right1 , name='Face_Right1')
    a.Set(faces=F_Left2 , name='Face_Left2')
    a.Set(faces=F_Right2 , name='Face_Right2')
    
    return ;


def Faces_PBC_1(elemsize):
    
    a = mdb.models['Model-%d' %(q)].rootAssembly     
    
    # Left-Right faces
    fLeftnodes=a.sets['Face_Left1'].nodes
    fRightnodes=a.sets['Face_Right1'].nodes
    
    fLeftCoord=[]
    fRightCoord=[]
    
    for node in fLeftnodes:
        fLeftCoord=fLeftCoord+[[node.coordinates[0],node.coordinates[1], node.coordinates[2],node.label]]        
    for node in fRightnodes:
        fRightCoord=fRightCoord+[[node.coordinates[0],node.coordinates[1], node.coordinates[2],node.label]]
        
    fLeftCoord.sort()
    fRightCoord.sort()
    
    NumLeft=len(fLeftCoord)
    p=a.instances['Part-1-1']
    Node_Tol=1e-10
    for i in range(0,NumLeft):
        if abs(fLeftCoord[i][0]-fRightCoord[i][0])<Node_Tol and abs(fLeftCoord[i][1]-fRightCoord[i][1])<Node_Tol:
            NLabel=fLeftCoord[i][3]
            a.Set(nodes=p.nodes[NLabel-1:NLabel], name='Face_Left1_'+str(i)) 
            NLabel=fRightCoord[i][3]
            a.Set(nodes=p.nodes[NLabel-1:NLabel], name='Face_Right1_'+str(i))    
        else:
            print('Distance between nodes is more than tolerance')
            break        

    for i in range(0,NumLeft):
        mdb.models['Model-%d' %(q)].Equation(name='Const-Faces_RightLeft1-x'+str(i), 
            terms=((1,'Face_Left1_'+str(i),1),(-1,'Face_Right1_'+str(i),1)))                    
        
        mdb.models['Model-%d' %(q)].Equation(name='Const-Faces_RightLeft1-y'+str(i), 
            terms=((1,'Face_Left1_'+str(i),2),(-1,'Face_Right1_'+str(i),2)))
        
        mdb.models['Model-%d' %(q)].Equation(name='Const-Faces_RightLeft1-z'+str(i), 
            terms=((1,'Face_Left1_'+str(i),3),(-1,'Face_Right1_'+str(i),3)))     
        
    return ;



def Faces_PBC_2(elemsize, W_sample, W_connecting1, W_inf):
    
    a = mdb.models['Model-%d' %(q)].rootAssembly     
    
    # Left-Right faces
    fLeftnodes=a.sets['Face_Left2'].nodes
    fRightnodes=a.sets['Face_Right2'].nodes
    
    fLeftCoord=[]
    fRightCoord=[]
    
    for node in fLeftnodes:
        fLeftCoord=fLeftCoord+[[node.coordinates[0],node.coordinates[1], node.coordinates[2],node.label]]        
    for node in fRightnodes:
        fRightCoord=fRightCoord+[[node.coordinates[0],node.coordinates[1], node.coordinates[2],node.label]]
        
    fLeftCoord.sort()
    fRightCoord.sort()
    
    NumLeft=len(fLeftCoord)
    p=a.instances['Part-1-1']
    Node_Tol=1e-10
    for i in range(0,NumLeft):
        
        if fLeftCoord[i][0] < W_sample + W_connecting1 + W_inf / 2:
            continue 
                
        if abs(fLeftCoord[i][0]-fRightCoord[i][0])<Node_Tol and abs(fLeftCoord[i][1]-fRightCoord[i][1])<Node_Tol:
            NLabel=fLeftCoord[i][3]
            a.Set(nodes=p.nodes[NLabel-1:NLabel], name='Face_Left2_'+str(i)) 
            NLabel=fRightCoord[i][3]
            a.Set(nodes=p.nodes[NLabel-1:NLabel], name='Face_Right2_'+str(i))    
        else:
            print('Distance between nodes is more than tolerance')
            break        

    for i in range(0,NumLeft):
        
        if fLeftCoord[i][0] < W_sample + W_connecting1 + W_inf / 2:
            continue 
                          
        mdb.models['Model-%d' %(q)].Equation(name='Const-Faces_RightLeft2-x'+str(i), 
            terms=((1,'Face_Left2_'+str(i),1),(-1,'Face_Right2_'+str(i),1)))                    
        
        mdb.models['Model-%d' %(q)].Equation(name='Const-Faces_RightLeft2-y'+str(i), 
            terms=((1,'Face_Left2_'+str(i),2),(-1,'Face_Right2_'+str(i),2)))
        
        mdb.models['Model-%d' %(q)].Equation(name='Const-Faces_RightLeft2-z'+str(i), 
            terms=((1,'Face_Left2_'+str(i),3),(-1,'Face_Right2_'+str(i),3)))     
        
    return ;
###############################################################################

##----------------------Set the background in ABAQUS white---------------------
###############################################################################
session.graphicsOptions.setValues(backgroundStyle=SOLID, 
    backgroundColor='#FFFFFF', translucencyMode=2)

session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    datumAxes=OFF, datumPlanes=OFF)
    
session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
    datumAxes=OFF, datumPlanes=OFF)
###############################################################################

##----------------------------Pulse Parameters---------------------------------
###############################################################################
f=5e6
Nc=5
wf=2*np.pi*f
T=1/f
Tc=Nc*T
t=np.linspace(0,Tc,10000)
t=t.reshape(t.shape[0],1)
time=t
g1=np.sin(wf*time)

start_time = 5e-9
end_time = 1.6e-5#1.95e-5#0.92e-5#1.9e-5

time = np.append(time, end_time) + start_time
g1 = np.append(g1, 0.0)

time = np.append(0.0, time)
g1 = np.append(0.0, g1)
l1 = g1

time=time.reshape(time.shape[0],1)
l1=l1.reshape(l1.shape[0],1)
###############################################################################

##-----------------------Material Properties for sample------------------------
###############################################################################

E_sample=72.997166e9
v_sample=0.335345
dens_sample=2579.433958406887

alpha_sample=0.0
beta_sample=0.4e-11 # 1e-11, 8e-12, 4e-12
cp_sample = np.sqrt(E_sample*(1-v_sample)/((1+v_sample)*(1-2*v_sample)*dens_sample))

wavelength = cp_sample / f
###############################################################################

##----------------------------Sample Dimensions--------------------------------
###############################################################################
H_sample = 10.19e-3#10.217e-3
W_sample = 0.0016526605425717712#1.664e-3/1.0 # 0.0012889288576178283, 0.0012889288576178283 * 1.6666666666666667
L_sample = 0.00163#1.664e-3/1.0
###############################################################################

##----------------------------Porosity Parameters------------------------------
###############################################################################
porosity_range = 4 * np.linspace(4, 10, 1)#[0, 38, 76, 114]
#rmean_range=list(np.linspace(30e-6,100e-6,int((100e-6-30e-6)/5e-6+1.0)))
#std_range=list(np.linspace(15e-6,25e-6, 3))

rmean_range =[220e-6, ] # 100e-6 #213.67819287444507e-6
std_range=[0.0,]
###############################################################################

##-----------------Element and Infinite Layer size calculation-----------------
###############################################################################
elemsize_range=[wavelength / 40, ]#[wavelength / 10, wavelength / 15, wavelength / 20, wavelength / 25, 
#                wavelength / 30, wavelength / 35, wavelength / 40, wavelength / 45,
 #               wavelength / 50, wavelength / 55, wavelength / 60, wavelength / 65]# 40 is good for HEX and 60 for TET (maybe try 70 for HEX also)
elemType4 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT,
                          secondOrderAccuracy=OFF, distortionControl=DEFAULT)               
elemType8 = mesh.ElemType(elemCode=AC3D8, secondOrderAccuracy=OFF, elemLibrary=EXPLICIT)


wavelength_div_range = [1/1.5, ]
###############################################################################

q = 6
for elemsize in elemsize_range:
    for wavelength_div in wavelength_div_range:
        N_Hex = 3#int(0.35e-3 / elemsize)
        
    #    if elemsize > 2 * np.min(rmean_range):
    #        elemsize = np.min(rmean_range) / 2
        
        W_connecting = N_Hex * elemsize
        H_connecting = N_Hex * elemsize
        L_connecting = N_Hex * elemsize
        W_connecting1 = wavelength / wavelength_div
        W_inf = 30 * elemsize
        for porosity in porosity_range:
            for rmean in rmean_range:
                for std in std_range:
                    
                    q=q+1
                    k=0
                    while k < 1:
                        k = k + 1
                        
                        populatedvolume=0
##-----------------Calculation of coordinates and radius of pores-------------
###############################################################################
                        if int(porosity)>0.0:
                            num_pores=0
                            rmin=rmean-std
        
                            ic=0
                            
                            ri1=np.array([])
                            xi1=np.array([])
                            yi1=np.array([])
                            zi1=np.array([])
                            
                            cr=np.ones((int(porosity))) * rmean
                            
                            i=0
                            j=0
                            while i<10000 and num_pores < porosity:
                                i=i+1
                                xmin=rmin
                                ymin=rmin
                                zmin=rmin
                                xmax=W_sample+W_connecting1-rmin
                                ymax=H_sample-rmin
                                zmax=L_sample-rmin
                            
                                r = cr[j]
                                x=xmin+(xmax-xmin)*random.random()
                                y=ymin+(ymax-ymin)*random.random()
                                z=zmin+(zmax-zmin)*random.random()
                            
                                A=CheckBorder(x,y,z,r,elemsize, W_connecting, H_connecting, L_connecting,W_sample+W_connecting1,H_sample,L_sample)
                            
                                if A==True:
                                    if ic==0:
                                        ic=1
                                        B=True
                                    else:
                                        B=CheckOverlap(x,y,z,r,elemsize,xi1,yi1,zi1,ri1)
                                        if B==True:
                                            ic=ic+1
                                            if A and B == True:
                                                num_pores=num_pores+1
                                                populatedvolume=populatedvolume+4*np.pi*(r**3)/3.0
                                                ri1=np.append(ri1,r)
                                                xi1=np.append(xi1,x)
                                                yi1=np.append(yi1,y)
                                                zi1=np.append(zi1,z)
                                                j=j+1
###############################################################################

##------------------------------Creation of Parts------------------------------
###############################################################################
                        mdb.Model(name='Model-%d' %(q), modelType=STANDARD_EXPLICIT)
                
                        mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', sheetSize=W_sample)        
                        mdb.models['Model-%d' %(q)].sketches['__profile__'].rectangle(point1=(0,0),point2=(W_sample + W_connecting1 + W_inf,H_sample))
                        
                        mdb.models['Model-%d' %(q)]. Part ( dimensionality = THREE_D , name ='Part-1'
                            , type = 
                            DEFORMABLE_BODY)
                        mdb.models['Model-%d' %(q)].parts ['Part-1'].BaseSolidExtrude(sketch =
                            mdb.models['Model-%d' %(q)].sketches['__profile__'], depth = L_sample)
                        del mdb.models['Model-%d' %(q)].sketches['__profile__']                                 
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 2
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=W_sample + W_connecting1)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 4
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=L_connecting)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 6
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=L_sample - L_connecting)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 8
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=H_connecting)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 10
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=H_sample - H_connecting)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 12
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=W_sample + W_connecting1 - W_connecting)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##----------------Creation of partition for infinite layer in z-axis-----------
###############################################################################
                        k = 14
                            
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=W_connecting)
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        c = p.cells
                        pickedCells = c.getSequenceFromMask(mask=('[#2 ]', ), )
                        d2 = p.datums
                        p.PartitionCellByDatumPlane(datumPlane=d2[k], cells=pickedCells)
###############################################################################

##------------------------------Creation of void-------------------------------
###############################################################################
                        if porosity>0.0:
                            for i in range(len(xi1)):
                                if i<492:
                                    k = 14 + 2*(i+1)
                                else:
                                    k = 14 + 2*(i+1) + 1
                                mdb.models['Model-%d' %(q)].parts['Part-1'].DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=zi1[i])
        
                                f1, e1, d2 = mdb.models['Model-%d' %(q)].parts ['Part-1'].faces, mdb.models['Model-%d' %(q)].parts ['Part-1'].edges, mdb.models['Model-%d' %(q)].parts ['Part-1'].datums
                                t = mdb.models['Model-%d' %(q)].parts ['Part-1'].MakeSketchTransform(sketchPlane=d2[k], sketchUpEdge=e1[63], 
                                        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, zi1[i]))
                                
                                mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', sheetSize=W_sample, transform=t)
                                mdb.models['Model-%d' %(q)].sketches['__profile__'].ArcByCenterEnds(center=(xi1[i], yi1[i]), point1=(xi1[i], yi1[i]+ri1[i]), point2=(xi1[i], yi1[i]-ri1[i]), direction=CLOCKWISE)
                                mdb.models['Model-%d' %(q)].sketches['__profile__'].Line(point1=(xi1[i], yi1[i]+ri1[i]), point2=(xi1[i], yi1[i]-ri1[i]))
                                mdb.models['Model-%d' %(q)].sketches['__profile__'].ConstructionLine(point1=(xi1[i], yi1[i]+ri1[i]), point2=(xi1[i], yi1[i]-ri1[i]))
                   
                                f, e, d1 = mdb.models['Model-%d' %(q)].parts ['Part-1'].faces, mdb.models['Model-%d' %(q)].parts ['Part-1'].edges, mdb.models['Model-%d' %(q)].parts ['Part-1'].datums
                    
                                mdb.models['Model-%d' %(q)].parts ['Part-1'].CutRevolve(sketchPlane=d1[k], sketchUpEdge=e[63], sketchPlaneSide=SIDE1, 
                                                                                        sketchOrientation=RIGHT, sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'], 
                                                                                        angle=360.0, flipRevolveDirection=OFF)
                                del mdb.models['Model-%d' %(q)].sketches['__profile__']
###############################################################################

##------------------------------Partition for load-----------------------------
###############################################################################
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        f1, e, d = p.faces, p.edges, p.datums
                        t = p.MakeSketchTransform(sketchPlane=f1[23], sketchUpEdge=e[60], 
                            sketchPlaneSide=SIDE1, origin=(0.0, H_sample, 0.0))
                        s = mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', 
                            sheetSize=W_sample, transform=t)
                        s.rectangle(point1=(0.0, 0.0), point2=(W_sample, -L_sample))
                    
                        p = mdb.models['Model-%d' %(q)].parts['Part-1']
                        f = p.faces
                        pickedFaces = f.getSequenceFromMask(mask=('[#4800000 #10 ]', ), )
                        e1, d2 = p.edges, p.datums
                        p.PartitionFaceBySketch(sketchUpEdge=e1[60], faces=pickedFaces, sketch=s)
                        del mdb.models['Model-%d' %(q)].sketches['__profile__']
###############################################################################

#-----------------Creation and assignment of material for sample---------------
###############################################################################
                        mdb.models['Model-%d' %(q)].Material(name='Material_sample')
                        mdb.models['Model-%d' %(q)].materials['Material_sample'].Elastic(table=((E_sample, v_sample), ))
                        mdb.models['Model-%d' %(q)].materials['Material_sample'].Density(table=((dens_sample, ), ))
                        mdb.models['Model-%d' %(q)].materials['Material_sample'].Damping(alpha_sample, beta_sample)
                    
                        mdb.models['Model-%d' %(q)].HomogeneousSolidSection(name='Section_sample',
                            material='Material_sample', thickness=None)
                        
                        cells = mdb.models['Model-%d' %(q)].parts['Part-1'].cells
                        region = (cells,)
                        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(region=region, sectionName='Section_sample', offset=0.0,
                            offsetType=MIDDLE_SURFACE, offsetField='',
                            thicknessAssignment=FROM_SECTION)                  
###############################################################################

#----------------------------Assembly creation---------------------------------
############################################################################### 
                        mdb.models['Model-%d' %(q)].rootAssembly.Instance(name='Part-1-1', part=mdb.models['Model-%d' %(q)].parts['Part-1'], dependent=OFF)
                                                
                        # The view inside the model is set to assembly
                        a = mdb.models['Model-%d' %(q)].rootAssembly
                        session.viewports['Viewport: 1'].setValues(displayedObject=a)
###############################################################################

#------------------------Detection of faces for PBC----------------------------
###############################################################################                
                        Faces_detection(H_sample, W_sample, W_connecting1, W_inf, L_sample)
###############################################################################

#------------Detection of edges of sample for infinite layer meshing-----------
###############################################################################     
                        a=mdb.models['Model-%d' %(q)].rootAssembly
                        
                        e=a.instances['Part-1-1'].edges 
                        
                        Edge_z1=e.findAt((W_sample + W_connecting1 + W_inf / 2, 0.0, 0.0))
                        Edge_z2=e.findAt((W_sample + W_connecting1 + W_inf / 2, 0.0, L_sample))
                        Edge_z3=e.findAt((W_sample + W_connecting1 + W_inf / 2, H_sample, 0.0))
                        Edge_z4=e.findAt((W_sample + W_connecting1 + W_inf / 2, H_sample, L_sample))
                        
                        q1=Edge_z1.index
                        q2=Edge_z2.index 
                        q3=Edge_z3.index 
                        q4=Edge_z4.index 
                        
                        E_z1=e[q1:q1+1]
                        E_z2=e[q2:q2+1]
                        E_z3=e[q3:q3+1]
                        E_z4=e[q4:q4+1]
                
                        a.Set(edges=E_z1 , name='Ez1')
                        a.Set(edges=E_z2 , name='Ez2')
                        a.Set(edges=E_z3 , name='Ez3')
                        a.Set(edges=E_z4 , name='Ez4')
###############################################################################

#-----------------------------Seeding for meshing------------------------------
###############################################################################
                        a=mdb.models['Model-%d' %(q)].rootAssembly
                                            
                        a.seedEdgeByNumber(edges=E_z1, number=1, constraint=FIXED)
                        a.seedEdgeByNumber(edges=E_z2, number=1, constraint=FIXED)
                        a.seedEdgeByNumber(edges=E_z3, number=1, constraint=FIXED)
                        a.seedEdgeByNumber(edges=E_z4, number=1, constraint=FIXED)
                        
                        a=mdb.models['Model-%d' %(q)].rootAssembly  
                        
                        partInstances=(a.instances['Part-1-1'],)
                        a.seedPartInstance(regions=partInstances, size=elemsize, deviationFactor=0.1, minSizeFactor=0.1)
###############################################################################

#-----------------------------Assignment of elements types---------------------
###############################################################################
                        order = 1
                        if order == 1:
                            elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT,
                                                      secondOrderAccuracy=OFF, distortionControl=DEFAULT)
                            elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
                            elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT,
                                                      secondOrderAccuracy=OFF, distortionControl=DEFAULT)
                           
                        
                        a=mdb.models['Model-%d' %(q)].rootAssembly    
        
                        c1 = a.instances['Part-1-1'].cells.findAt(((W_connecting + elemsize / 4, H_connecting + elemsize / 4, L_connecting + elemsize / 4),)) 
                        if porosity>0.0:
                            pickedRegions=c1
                            mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
                            pickedRegions=(c1,)
                            a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3,)) 
                        elif porosity == 0.0:
                            pickedRegions=c1
                            mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
                            pickedRegions=(c1, )
                            a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3,))

#------------------------------------------------------------------------------
                        c2 = a.instances['Part-1-1'].cells.findAt(((W_sample + W_connecting1 + W_inf, 0.0, 0.0),))
                        pickedRegions=c2
                        a.setMeshControls(regions=pickedRegions, technique=SWEEP, 
                            algorithm=ADVANCING_FRONT)
                        pickedRegions=(c2,)
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType8,))
        
                        c2 = a.instances['Part-1-1'].cells
                        pickedCells = c2.getSequenceFromMask(mask=('[#80 ]', ), )
                        e1 = a.instances['Part-1-1'].edges
                        f1 = a.instances['Part-1-1'].faces
                        a.assignStackDirection(referenceRegion=f1[37], cells=pickedCells)
#------------------------------------------------------------------------------
    
                        c3 = a.instances['Part-1-1'].cells.findAt(((W_sample + W_connecting1 - elemsize / 4, H_sample / 2,  L_sample / 2),))
                        pickedRegions=c3
                        mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=HEX, technique=STRUCTURED)
                        pickedRegions=(c3, )
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType4, ))
    
                        c4 = a.instances['Part-1-1'].cells.findAt(((elemsize / 4, H_sample / 2,  L_sample / 2),))
                        pickedRegions=c4
                        mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=HEX, technique=STRUCTURED)
                        pickedRegions=(c4, )
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType4, ))
                        
                        c5 = a.instances['Part-1-1'].cells.findAt(((W_sample / 2, H_sample / 2,  elemsize / 4),))
                        pickedRegions=c5
                        mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=HEX, technique=STRUCTURED)
                        pickedRegions=(c5, )
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType4, ))
    
                        c6 = a.instances['Part-1-1'].cells.findAt(((W_sample / 2, H_sample / 2,  L_sample - elemsize / 4),))
                        pickedRegions=c6
                        mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=HEX, technique=STRUCTURED)
                        pickedRegions=(c6, )
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType4, ))   
                        
                        c7 = a.instances['Part-1-1'].cells.findAt(((W_sample / 2, elemsize / 4,  L_sample / 2),))
                        pickedRegions=c7
                        mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=HEX, technique=STRUCTURED)
                        pickedRegions=(c7, )
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType4, ))      
                        
                        c8 = a.instances['Part-1-1'].cells.findAt(((W_sample / 2, H_sample -  elemsize / 4,  L_sample / 2),))
                        pickedRegions=c8
                        mdb.models['Model-%d' %(q)].rootAssembly.setMeshControls(regions=pickedRegions, elemShape=HEX, technique=STRUCTURED)
                        pickedRegions=(c8, )
                        a.setElementType(regions=pickedRegions, elemTypes=(elemType4, ))    
        
                        partInstances =(a.instances['Part-1-1'], )
                        a.generateMesh(regions=partInstances)
                        quality=a.verifyMeshQuality(ANALYSIS_CHECKS)
                        
                        if len(mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].nodes) != 0:
                            if len(quality["failedElements"])==0:#and len(quality["warningElements"])==0:                            
                                break
###############################################################################

#--------------------------Tabular Amplitude creation--------------------------
###############################################################################    
                    mdb.models['Model-%d' %(q)].TabularAmplitude(name='Amp-1', timeSpan=STEP,
                        smooth=SOLVER_DEFAULT, data=tuple(map(tuple, np.concatenate((time,l1),axis=1))))
###############################################################################

#---------------------------------Steps creation-------------------------------
###############################################################################
                    tcrit = elemsize / cp_sample
                    tincr = tcrit/4# tcrit / 3, tcrit/5
                    
                    tstep = Tc + 2 * start_time
                    mdb.models['Model-%d' %(q)].ExplicitDynamicsStep(name='Step-1', previous='Initial',
                        timePeriod=tstep, timeIncrementationMethod=FIXED_USER_DEFINED_INC,
                        userDefinedInc=tincr, linearBulkViscosity=0.06, quadBulkViscosity=1.2, improvedDtMethod=ON)
    
                    tstep = end_time-tstep
                    mdb.models['Model-%d' %(q)].ExplicitDynamicsStep(name='Step-2', previous='Step-1',
                        timePeriod=tstep, timeIncrementationMethod=FIXED_USER_DEFINED_INC,
                        userDefinedInc=tincr, linearBulkViscosity=0.06, quadBulkViscosity=1.2, improvedDtMethod=ON)
###############################################################################

#---------------------Creation of symmetry condition in x-axis - 1-------------
###############################################################################                
                    xMin= - elemsize/4.0 - 1e-10
                    xMax= elemsize/4.0 + 1e-10
                    
                    yMin= - elemsize/4.0 - 1e-10
                    yMax= H_sample + elemsize/4.0 +1e-10
                    
                    zMin= - elemsize/4.0 - 1e-10
                    zMax= L_sample + elemsize/4.0 + 1e-10
                    
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.getByBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)
    
                    regionDef= regionToolset.Region(faces=v1)
                    mdb.models['Model-%d' %(q)].XsymmBC(name='XSYMM-BC-1', createStepName='Initial', 
                        region=regionDef, localCsys=None)
###############################################################################

#---------------------Creation of fixture in y-axis - 1------------------------
################################################################################                
                    xMin= - elemsize/4.0 - 1e-10
                    xMax= elemsize/4.0 + 1e-10
                    
                    yMin= - elemsize/4.0 - 1e-10
                    yMax= elemsize/4.0 +1e-10
                    
                    zMin= - elemsize/4.0 - 1e-10
                    zMax= elemsize/4.0 + 1e-10
                    
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].vertices.getByBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)
    
                    regionDef= regionToolset.Region(vertices=v1)
                    mdb.models['Model-%d' %(q)].rootAssembly.Set(name='FixtureNode-1', region=regionDef)
                    
                    regionDef = mdb.models['Model-%d' %(q)].rootAssembly.sets['FixtureNode-1']
                    mdb.models['Model-%d' %(q)].DisplacementBC(name='Fixture-1', createStepName='Initial', 
                        region=regionDef, u1=UNSET, u2=SET, u3=UNSET, amplitude=UNSET, 
                        distributionType=UNIFORM, fieldName='', localCsys=None)
################################################################################

#---------------------Creation of fixture in y-axis - 2------------------------
################################################################################                
                    xMin= W_sample + W_connecting1 + W_inf - elemsize/4.0 - 1e-10
                    xMax= W_sample + W_connecting1 + W_inf+ elemsize/4.0 + 1e-10
                    
                    yMin= - elemsize/4.0 - 1e-10
                    yMax= elemsize/4.0 +1e-10
                    
                    zMin= - elemsize/4.0 - 1e-10
                    zMax= elemsize/4.0 + 1e-10
                    
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].vertices.getByBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)
    
                    regionDef= regionToolset.Region(vertices=v1)
                    mdb.models['Model-%d' %(q)].rootAssembly.Set(name='FixtureNode-2', region=regionDef)
                    
                    regionDef = mdb.models['Model-%d' %(q)].rootAssembly.sets['FixtureNode-2']
                    mdb.models['Model-%d' %(q)].DisplacementBC(name='Fixture-2', createStepName='Initial', 
                        region=regionDef, u1=UNSET, u2=SET, u3=UNSET, amplitude=UNSET, 
                        distributionType=UNIFORM, fieldName='', localCsys=None)
################################################################################

#---------------------Creation of fixture in y-axis - 3------------------------
################################################################################                
                    xMin= - elemsize/4.0 - 1e-10
                    xMax= elemsize/4.0 + 1e-10
                    
                    yMin= - elemsize/4.0 - 1e-10
                    yMax= elemsize/4.0 +1e-10
                    
                    zMin= L_sample - elemsize/4.0 - 1e-10
                    zMax= L_sample + elemsize/4.0 + 1e-10
                    
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].vertices.getByBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)
    
                    regionDef= regionToolset.Region(vertices=v1)
                    mdb.models['Model-%d' %(q)].rootAssembly.Set(name='FixtureNode-3', region=regionDef)
                    
                    regionDef = mdb.models['Model-%d' %(q)].rootAssembly.sets['FixtureNode-3']
                    mdb.models['Model-%d' %(q)].DisplacementBC(name='Fixture-3', createStepName='Initial', 
                        region=regionDef, u1=UNSET, u2=SET, u3=UNSET, amplitude=UNSET, 
                        distributionType=UNIFORM, fieldName='', localCsys=None)
################################################################################

#---------------------Creation of fixture in y-axis - 4------------------------
################################################################################                
                    xMin= W_sample + W_connecting1 + W_inf - elemsize/4.0 - 1e-10
                    xMax= W_sample + W_connecting1 + W_inf+ elemsize/4.0 + 1e-10
                    
                    yMin= - elemsize/4.0 - 1e-10
                    yMax= elemsize/4.0 +1e-10
                    
                    zMin= L_sample - elemsize/4.0 - 1e-10
                    zMax= L_sample + elemsize/4.0 + 1e-10
                    
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].vertices.getByBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)
    
                    regionDef= regionToolset.Region(vertices=v1)
                    mdb.models['Model-%d' %(q)].rootAssembly.Set(name='FixtureNode-4', region=regionDef)
                    
                    regionDef = mdb.models['Model-%d' %(q)].rootAssembly.sets['FixtureNode-4']
                    mdb.models['Model-%d' %(q)].DisplacementBC(name='Fixture-4', createStepName='Initial', 
                        region=regionDef, u1=UNSET, u2=SET, u3=UNSET, amplitude=UNSET, 
                        distributionType=UNIFORM, fieldName='', localCsys=None)
################################################################################

#------------------------------Creation of PBC---------------------------------
###############################################################################
                    Faces_PBC_1(elemsize)
                    Faces_PBC_2(elemsize, W_sample, W_connecting1, W_inf)
###############################################################################

#------------------------Creation of displacement load-------------------------
###############################################################################
                    xMin= - elemsize/4.0 - 1e-10
                    xMax= W_sample + elemsize/4.0 + 1e-10
                    
                    yMin= H_sample - elemsize/4.0 - 1e-10
                    yMax= H_sample + elemsize/4.0 +1e-10
                    
                    zMin= - elemsize/4.0 - 1e-10
                    zMax= L_sample + elemsize/4.0 + 1e-10
                
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.getByBoundingBox(xMin,yMin,zMin,xMax,yMax,zMax)
    
                    regionDef= regionToolset.Region(faces=v1)  
                    
                    mdb.models['Model-%d' %(q)].DisplacementBC(name='Load-1', createStepName='Step-1', 
                        region=regionDef, u1=UNSET, u2=-1.9e-15, u3=UNSET, amplitude='Amp-1', 
                        fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
                    mdb.models['Model-%d' %(q)].boundaryConditions['Load-1'].deactivate('Step-2') 
###############################################################################

#-------------------------Set_creation to record signal------------------------
###############################################################################
                    xMin= - elemsize/4.0 - 1e-10
                    xMax= W_sample + elemsize/4.0 + 1e-10
                    
                    yMin= H_sample - elemsize/4.0 - 1e-10
                    yMax= H_sample + elemsize/4.0 +1e-10
                    
                    zMin= - elemsize/4.0 - 1e-10
                    zMax= L_sample + elemsize/4.0 + 1e-10
                    
                    v1=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.getByBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)
                    
                    regionDef= regionToolset.Region(faces=v1)
                    mdb.models['Model-%d' %(q)].rootAssembly.Set(name='RE1', region=regionDef) 
###############################################################################

#-------------------Field output creation to record signal---------------------
###############################################################################
                    regionDef=mdb.models['Model-%d' %(q)].rootAssembly.sets['RE1']
                    frequency = int(np.ceil((1/250e6) / tincr))
                    mdb.models['Model-%d' %(q)].HistoryOutputRequest(name='H-Output-2', 
                        createStepName='Step-1', variables=('U2',), frequency=frequency, 
                        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)  
###############################################################################   

#------------------------------Job creation------------------------------------
###############################################################################
                    mdb.Job(name='Job-%d' %(q), model='Model-%d' %(q), description='', type=ANALYSIS, 
                        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK, 
                        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, 
                        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=1, 
                        activateLoadBalancing=True, numThreadsPerMpiProcess=1, 
                        multiprocessingMode=DEFAULT, numCpus=1)  
                    
                    mdb.jobs['Job-%d' %(q)].writeInput(consistencyChecking=OFF)
###############################################################################

#----------Replacement of Acoustic with Infinite elements in the .inp file-----
###############################################################################
                    with open('Job-%d.inp' %(q), 'r') as file :
                      filedata = file.read()
                
                    filedata = filedata.replace('AC3D8', 'CIN3D8')
                    
                    with open('Job-%d.inp' %(q), 'w') as file:
                      file.write(filedata)
###############################################################################




