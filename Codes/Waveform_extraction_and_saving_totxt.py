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


#------------------------------------------------------------------------------
qq = [1, ]#[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
for q in qq:

    INPfilename='Job-%d' %(q)       # set odb name here
    path='./'                       # set odb path here (if in working dir no need to change!)
    INPpath=path+INPfilename+'.inp'      

    mdb.ModelFromInputFile(inputFileName=INPpath
        , name='Model-%d' %(q) )

    a=mdb.models['Model-%d' %(q)].rootAssembly    
    REnodes=a.sets['RE1'].nodes
    
    path='./'  
    file_name='Job-%d' %(q)
    
    o1 = session.openOdb('{}{}.odb'.format(path, file_name))
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    odb = session.odbs['{}{}.odb'.format(path, file_name)]
    
    data=0.0
    i=0.0
    for node in REnodes:
        i=i+1.0
        history_output='Spatial displacement: U2 at Node {} in NSET RE1'.format(node.label)
        xy = xyPlot.XYDataFromHistory(odb=odb, outputVariableName=history_output, steps=('Step-1','Step-2'), suppressQuery=True, __linkedVpName__='Viewport: 1')
        data=data+np.array(xy)
    
    data=data/i
    output_file_name='{}result_{}_U2.txt'.format(path,file_name)
    np.savetxt(output_file_name, data)

       
    

