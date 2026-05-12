# Visualize Atoms

これはase.Atomsの可視化のためのtips集です. **nglviewの出力はこのnotebookのpreviewには表示されないので注意してください.**


```python
import numpy as np
import nglview as nv
from ase.build import diamond100, molecule, bulk
from ase.io import read
from ase import Atoms
```


    



```python
def _get_standard_pos(atoms: Atoms) -&gt; np.ndarray:
    # NGLViewer shows atoms in standard positions.
    pos = atoms.positions
    if atoms.get_pbc().any():
        rcell, rot_t = atoms.cell.standard_form()
        standard_pos = pos.dot(rot_t.T)
    else:
        standard_pos = pos
    return standard_pos
```


## 分子と各原子のindexを同時に描画する関数
各原子の位置にその元素名とindexを描画する関数です.


```python
def view_with_atomindex(atoms: Atoms, label_color: str = "blue", label_scale: float = 1.0):
    v = nv.show_ase(atoms, viewer="ngl")
    v.add_label(
        color=label_color, labelType="text",
        labelText=[atoms[i].symbol + str(i) for i in range(atoms.get_global_number_of_atoms())],
        zOffset=1.0, attachment='middle_center', radius=label_scale
    )
    return v
```


```python
formula = "C60"
atoms = molecule(formula)
view_with_atomindex(atoms)
```


    NGLWidget()



```python
atoms = diamond100('C', size=(4,4,4),vacuum=10.0)
view_with_atomindex(atoms)
```


    NGLWidget()



```python
atoms = bulk('Pt', "fcc") * (4,4,4)
view_with_atomindex(atoms)
```


    NGLWidget()



## 各原子の座標を表示する関数

各原子の位置にマウスカーソルを合わせると座標が表示される関数です.


```python
darkgrey = [0.6, 0.6, 0.6]
def view_with_atomcoordinate(atoms: Atoms, radius=0.5):
    v = nv.show_ase(atoms, viewer="ngl")
    pos = atoms.positions
    standard_pos = _get_standard_pos(atoms)
    for i in range(len(atoms)):
        v.shape.add_sphere(standard_pos[i].tolist(), darkgrey, radius, f"x:{pos[i, 0]}, y:{pos[i, 1]}, z:{pos[i, 2]}")
    return v
```


```python
formula = "C60"
atoms = molecule(formula)
view_with_atomcoordinate(atoms)
```


    NGLWidget()



```python
atoms = diamond100('C', size=(4,4,4),vacuum=10.0)
view_with_atomcoordinate(atoms)
```


    NGLWidget()



```python
atoms = bulk('Pt', "fcc") * (4,4,4)
view_with_atomcoordinate(atoms)
```


    NGLWidget()



## 各原子の座標を表示し, indexと座標軸を描画する関数

`view_with_atomindex`, `view_with_atomcoordinate`の内容に加えてxyz座標軸を描画する関数です.


```python
darkgrey = [0.6, 0.6, 0.6]
grey = [0.8, 0.8, 0.8]
black = [0.0, 0.0, 0.0]

def view_with_axis(atoms, x_origin=None, y_origin=None, z_origin=None, length=None):
    pos = atoms.positions
    standard_pos = _get_standard_pos(atoms)

    min_pos = np.min(atoms.get_positions(), axis=0)
    max_diff = max(np.max(atoms.get_positions(), axis=0) - min_pos) + 0.5
    x_origin_auto, y_origin_auto, z_origin_auto =  min_pos - max_diff * 0.2
    length_auto = max_diff * 0.3
    x_origin = x_origin_auto if x_origin is None else x_origin
    y_origin = y_origin_auto if y_origin is None else y_origin
    z_origin = z_origin_auto if z_origin is None else z_origin
    length = length_auto if length is None else length

    v = nv.show_ase(atoms, viewer="ngl")
    v.shape.add_arrow([x_origin, y_origin, z_origin], [x_origin+length, y_origin, z_origin], grey, 0.3, "x-axis")
    v.shape.add_arrow([x_origin, y_origin, z_origin], [x_origin, y_origin+length, z_origin], grey, 0.3, "y-axis")
    v.shape.add_arrow([x_origin, y_origin, z_origin], [x_origin, y_origin, z_origin+length], grey, 0.3, "z-axis")
    v.shape.add_label([x_origin+length, y_origin, z_origin], black, 2, 'x')
    v.shape.add_label([x_origin, y_origin+length, z_origin], black, 2, 'y')
    v.shape.add_label([x_origin, y_origin, z_origin+length], black, 2, 'z')
    v.add_label(
        color="blue", labelType="text",
        labelText=[atoms[i].symbol + str(i) for i in range(atoms.get_global_number_of_atoms())],
        zOffset=1.0, attachment='middle_center', radius=1.0
    )
    for i in range(len(atoms)):
        v.shape.add_sphere(standard_pos[i].tolist(), darkgrey, 0.5, f"x:{pos[i, 0]}, y:{pos[i, 1]}, z:{pos[i, 2]}")
    return v
```


```python
formula = "C60"
atoms = molecule(formula)
view_with_axis(atoms)
```


    NGLWidget()



```python
atoms = diamond100('C', size=(4,4,4),vacuum=10.0)
view_with_axis(atoms)
```


    NGLWidget()



```python
atoms = bulk('Pt', "fcc") * (4,4,4)
view_with_axis(atoms)
```


    NGLWidget()

