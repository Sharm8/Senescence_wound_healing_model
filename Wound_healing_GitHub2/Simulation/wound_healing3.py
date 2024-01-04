
from cc3d import CompuCellSetup
        

from wound_healing3Steppables import MechanismsSteppable
CompuCellSetup.register_steppable(steppable=MechanismsSteppable(frequency=1))
       
        
from wound_healing3Steppables import PlotsSteppable
CompuCellSetup.register_steppable(steppable=PlotsSteppable(frequency=1))
     
  
from wound_healing3Steppables import ClearSteppable
CompuCellSetup.register_steppable(steppable=ClearSteppable(frequency=1))

      
from wound_healing3Steppables import WoundSteppable
CompuCellSetup.register_steppable(steppable=WoundSteppable(frequency=1))

        
from wound_healing3Steppables import MitosisSteppable
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))
     
        
from wound_healing3Steppables import SecretionSteppable
CompuCellSetup.register_steppable(steppable=SecretionSteppable(frequency=1))


#from wound_healing3Steppables import update_workspace_dirSteppable
#CompuCellSetup.register_steppable(steppable=update_workspace_dirSteppable(frequency=100))
       

CompuCellSetup.run()
