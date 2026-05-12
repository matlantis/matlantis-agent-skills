Copyright ENEOS, Corp., Preferred Computational Chemistry, Inc. and Preferred Networks, Inc. as contributors to Matlantis contrib project


###  原子数が多い構造のボンド表示対応（5000個程度まで） 

eneos ibuka 2022.12.15


```python
import ase,ase.io
from ase import Atoms 
import numpy as np
import nglview as nv

def add_lotsbonds(mols : ase.Atoms, view :  nv.widget.NGLWidget , rap_n = 200,  nglmax_n = 400 , lineOnly = False ):
    mols = mols[np.argsort(mols.positions[:,2])]
    myslices =  [mols[i : i+nglmax_n] for i in range(0, len(mols), nglmax_n-rap_n) ] 
    for myslice in myslices:
        structure = nv.ASEStructure(myslice)
        st = view.add_structure(structure,defaultRepresentation="")
        if lineOnly:
            st.add_ball_and_stick(lineOnly = True)
        else:
            st.add_ball_and_stick(cylinderOnly = True) #, linewidth=4

def lotsatoms_view(mols,radiusScale=0.3, rap_n = 200,  nglmax_n = 400 , lineOnly = False):
    """ NGLviewデフォルトでは400程度までしかbonds表示できない。分割して描画することで対応
        ｚ軸順に並べ替え、オーバラップさせながらボンド描画
    """    
    view = nv.NGLWidget(width=str(500)+ "px" ,height=str(500)+"px")
    structure = nv.ASEStructure(mols)
    #まずはspacefillのみ描く　
    view.add_structure(structure,defaultRepresentation="")
    view.add_spacefill(radiusType='covalent',radiusScale=radiusScale)
    view.center(component=0)
    view.add_unitcell()
    view.background ="#222222"
    view.camera  ="orthographic"
    add_lotsbonds(mols, view , rap_n = rap_n  , nglmax_n = nglmax_n , lineOnly = lineOnly )    
    return view
```


    



```python
# Ex1.　Cluster
from ase.cluster import wulff_construction

size = 2000  # Number of atoms
atoms = wulff_construction('Cu',  [(1, 0, 0), (1, 1, 0), (1, 1, 1)], [1.0, 1.1, 0.9] ,
                           size, 'fcc',rounding='above', latticeconstant=3.5)
atoms.numbers = ( list(range(24,30)) *1000 )[:len(atoms)]
print(len(atoms))
v = lotsatoms_view(atoms)
v
```

    2075



    NGLWidget(background='#222222')

