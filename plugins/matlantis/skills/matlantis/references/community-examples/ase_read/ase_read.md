Copyright Preferred Networks, Inc. as contributors to Matlantis contrib project

# ase.io.readがinitial_chargeを持つxyz fileのreadに失敗するbugに対するwork around
aseのversionによっては、initial_chargesがsetされたAtomsにchargesをsetしてase.io.write関数によってxyz fileに出力した際、ase.io.readでそのxyz fileをreadすることができないbugが確認されています。

このscriptでは上記bugの再現とそのwork aroundを紹介します。**なお、このbug fixは本家のase version3.23.0において適用され、解消されました。**

## bugが確認されているversionのaseをinstall


```python
# !pip install ase==3.22.1
```

## bugの再現
次のcellを実行するとinitial_chargesとchargesがsetされたatomsを保存したxyz fileをreadできないbugを再現できます


```python
import ase
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator


at = bulk('Si')
initial_charges = [1.0, -1.0]
charges = [-2.0, 2.0]

# initial_chargesをset
at.set_initial_charges(initial_charges)
at.calc = SinglePointCalculator(at, charges=charges)

# chargesをset
at.get_charges()

# xyz fileに出力
ase.io.write('charges1.xyz', at, format='extxyz')

# readに失敗する
r = ase.io.read('charges1.xyz')
```


    ---------------------------------------------------------------------------

    RuntimeError                              Traceback (most recent call last)

    &lt;ipython-input-2-6e19ffe14b16&gt; in &lt;module&gt;
         19 
         20 # readに失敗する
    ---&gt; 21 r = ase.io.read('charges1.xyz')
    

    ~/.local/lib/python3.7/site-packages/ase/io/formats.py in read(filename, index, format, parallel, do_not_split_by_at_sign, **kwargs)
        735     else:
        736         return next(_iread(filename, slice(index, None), format, io,
    --&gt; 737                            parallel=parallel, **kwargs))
        738 
        739 


    ~/.local/lib/python3.7/site-packages/ase/parallel.py in new_generator(*args, **kwargs)
        273             not kwargs.pop('parallel', True)):
        274             # Disable:
    --&gt; 275             for result in generator(*args, **kwargs):
        276                 yield result
        277             return


    ~/.local/lib/python3.7/site-packages/ase/io/formats.py in _iread(filename, index, format, io, parallel, full_output, **kwargs)
        801     # Make sure fd is closed in case loop doesn't finish:
        802     try:
    --&gt; 803         for dct in io.read(fd, *args, **kwargs):
        804             if not isinstance(dct, dict):
        805                 dct = {'atoms': dct}


    ~/.local/lib/python3.7/site-packages/ase/io/formats.py in wrap_read_function(read, filename, index, **kwargs)
        557         yield read(filename, **kwargs)
        558     else:
    --&gt; 559         for atoms in read(filename, index, **kwargs):
        560             yield atoms
        561 


    ~/.local/lib/python3.7/site-packages/ase/io/extxyz.py in read_xyz(fileobj, index, properties_parser)
        771         # check for consistency with frame index table
        772         assert int(fileobj.readline()) == natoms
    --&gt; 773         yield _read_xyz_frame(fileobj, natoms, properties_parser, nvec)
        774 
        775 


    ~/.local/lib/python3.7/site-packages/ase/io/extxyz.py in _read_xyz_frame(lines, natoms, properties_parser, nvec)
        507 
        508     for name, array in arrays.items():
    --&gt; 509         atoms.new_array(name, array)
        510 
        511     if duplicate_numbers is not None:


    ~/.local/lib/python3.7/site-packages/ase/atoms.py in new_array(self, name, a, dtype, shape)
        464 
        465         if name in self.arrays:
    --&gt; 466             raise RuntimeError('Array {} already present'.format(name))
        467 
        468         for b in self.arrays.values():


    RuntimeError: Array initial_charges already present


以下のように、initial chargesをexplicitにsetしない場合もreadの際にinitial chargesがsetされることでbugが再現します。


```python
import ase
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator


at = bulk('Si')
charges = [-2.0, 2.0]

at.calc = SinglePointCalculator(at, charges=charges)

# chargesをset
at.get_charges()

# xyz fileに出力
ase.io.write('charges2.xyz', at, format='extxyz')
# read
at = ase.io.read('charges2.xyz')
# chargesがinitial_chargesにsetされる
print("initial charges =", at.get_initial_charges())

at.calc = SinglePointCalculator(at, charges=charges)

# chargesをset
at.get_charges()

# xyz fileに出力
ase.io.write('charges2.xyz', at, format='extxyz')
# read
at = ase.io.read('charges2.xyz')
```

    initial charges = [-2.  2.]



    ---------------------------------------------------------------------------

    RuntimeError                              Traceback (most recent call last)

    &lt;ipython-input-3-5a8a61c64f9c&gt; in &lt;module&gt;
         27 ase.io.write('charges2.xyz', at, format='extxyz')
         28 # read
    ---&gt; 29 at = ase.io.read('charges2.xyz')
    

    ~/.local/lib/python3.7/site-packages/ase/io/formats.py in read(filename, index, format, parallel, do_not_split_by_at_sign, **kwargs)
        735     else:
        736         return next(_iread(filename, slice(index, None), format, io,
    --&gt; 737                            parallel=parallel, **kwargs))
        738 
        739 


    ~/.local/lib/python3.7/site-packages/ase/parallel.py in new_generator(*args, **kwargs)
        273             not kwargs.pop('parallel', True)):
        274             # Disable:
    --&gt; 275             for result in generator(*args, **kwargs):
        276                 yield result
        277             return


    ~/.local/lib/python3.7/site-packages/ase/io/formats.py in _iread(filename, index, format, io, parallel, full_output, **kwargs)
        801     # Make sure fd is closed in case loop doesn't finish:
        802     try:
    --&gt; 803         for dct in io.read(fd, *args, **kwargs):
        804             if not isinstance(dct, dict):
        805                 dct = {'atoms': dct}


    ~/.local/lib/python3.7/site-packages/ase/io/formats.py in wrap_read_function(read, filename, index, **kwargs)
        557         yield read(filename, **kwargs)
        558     else:
    --&gt; 559         for atoms in read(filename, index, **kwargs):
        560             yield atoms
        561 


    ~/.local/lib/python3.7/site-packages/ase/io/extxyz.py in read_xyz(fileobj, index, properties_parser)
        771         # check for consistency with frame index table
        772         assert int(fileobj.readline()) == natoms
    --&gt; 773         yield _read_xyz_frame(fileobj, natoms, properties_parser, nvec)
        774 
        775 


    ~/.local/lib/python3.7/site-packages/ase/io/extxyz.py in _read_xyz_frame(lines, natoms, properties_parser, nvec)
        507 
        508     for name, array in arrays.items():
    --&gt; 509         atoms.new_array(name, array)
        510 
        511     if duplicate_numbers is not None:


    ~/.local/lib/python3.7/site-packages/ase/atoms.py in new_array(self, name, a, dtype, shape)
        464 
        465         if name in self.arrays:
    --&gt; 466             raise RuntimeError('Array {} already present'.format(name))
        467 
        468         for b in self.arrays.values():


    RuntimeError: Array initial_charges already present


## xyz fileを保存する前のwork around
xyz fileを保存する前であれば、ase.io.write関数を実行する前にinitial_chargesをresetすることでreadできないbugを回避可能です。


```python
at = bulk('Si')
initial_charges = [1.0, -1.0]
charges = [-2.0, 2.0]

# initial_chargesをset
at.set_initial_charges(initial_charges)
at.calc = SinglePointCalculator(at, charges=charges)

# chargesをset
at.get_charges()

# initial_chargesをreset
at.set_initial_charges()

# xyz fileに出力
ase.io.write('charges3.xyz', at, format='extxyz')

# read
r = ase.io.read('charges3.xyz')

print(r)
# chargesがinitial_chargesにsetされる
print(r.get_initial_charges())
```

    Atoms(symbols='Si2', pbc=True, cell=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]], initial_charges=...)
    [-2.  2.]


## xyz file保存後にreadする方法
ase_ioにaseのread関数の実装に変更を加えたsource codeを配置しました。ase.io.readで読み込めなくなってしまったxyz fileも以下のread関数で読み込むことができます。


```python
from ase_io import read
```


```python
r = read('charges1.xyz')
# r = read('charges2.xyz')

print(r)
print(r.get_initial_charges())
print(r.get_charges())
```

    Atoms(symbols='Si2', pbc=True, cell=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]], charges=..., initial_charges=..., calculator=SinglePointCalculator(...))
    [ 1. -1.]
    [-2.  2.]

