import sys
from pyworkflow.project import *
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3 import *

path = sys.argv[1]

proj = Project(path)
proj.create()

from tests.tester import *

#prot = MyProtocol(name="test3", n=2, workingDir=proj.getPath('Runs', 'T3'))
#proj.launchProtocol(prot)

prot2 = ProtImportMicrographs(workingDir=proj.getPath('Runs', 'Import1'), 
                              pattern="InputData/*.mrc", sampling=1, voltage=200)
proj.launchProtocol(prot2)

l = proj.mapper.select(classname='SetOfMicrographs')
#l = proj.mapper.getAll()

for p in l:
    p.printAll()

prot3 = XmippProtDownsampleMicrographs(workingDir=proj.getPath('Runs', 'Downsample1'),
                                       inputMicrographs=l[0], downFactor=2)

proj.launchProtocol(prot3)

l = proj.mapper.select(classname='SetOfMicrographs')
#l = proj.mapper.getAll()

for p in l:
    p.printAll()
