Copyright ENEOS Corporation as contributors to Matlantis contrib project.

# Build genstr in ATAT package

This example depends on genstr in ATAT package.\
Please download ATAT from the home page, and put the tar.gz file in ./packages/ directory.

[ATAT HomePage](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/) -&gt; [Download:  Whole toolkit (Stable version)](http://alum.mit.edu/www/avdw/atat/atat3_36.tar.gz)


```python
import os
import glob
import shutil
import tarfile
from subprocess import run
from pathlib import Path
```

### Extract the package from the archive


```python
rootdir = Path().resolve()
print(rootdir)
```

    /home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example



```python
!ls -l packages/atat*.*
```

    -rw-r--r-- 1 jovyan 1000 6993920 Feb  7 00:25 packages/atat3_36.tar.gz



```python
destdir = 'packages'
os.makedirs(destdir, exist_ok=True)
tgzfile = glob.glob(f'{destdir}/atat*.tar.gz')[0]

tar = tarfile.open(tgzfile)
tar.extractall(path='packages/')
```


```python
atatroot = Path(glob.glob(f'{destdir}/atat')[0]).resolve()
atatroot
```




    PosixPath('/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/packages/atat')



### Build genstr


```python
os.chdir(rootdir)
shutil.copy(f'packages/patches/foolproof', f'{str(atatroot)}')
shutil.copy(f'packages/patches/makefile', f'{str(atatroot)}')
```




    '/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/packages/atat/makefile'




```python
os.chdir(f'{atatroot}')
cmd = f'make genstr'.split()
p = run(cmd, shell=False)
os.chdir(f'{rootdir}')
```

    bash ./foolproof /home/jovyan   "g++"
    I will now perform various checks to see if you have the appropriate software
    installed. To override these checks, type
    make force
    /home/jovyan g++
    Testing compiler features...
      : patched
      friend : patched
      doth : patched
      str2s : patched
      template : patched
      uT : ok
      static : patched
    Tests successfully passed.
    make -C src clean
    make[1]: Entering directory '/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/packages/atat/src'
    rm -f *.o *.cc *.c *.h *.exe
    rm -f maps mmaps emc2 phb checkcell corrdump kmesh genstr gensqs mcsqs nntouch fixcell csfit cv cellcvrt lsfit fitsvsl svsl felec pdef fitfc nnshell memc2 gce gencs triph strpath infdet apb
    make[1]: Leaving directory '/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/packages/atat/src'
    make -C src genstr   "CXX=g++"
    make[1]: Entering directory '/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/packages/atat/src'
    ./patchlang &lt; anyfft.hh &gt; anyfft.h
    ./patchlang &lt; drawpd.hh &gt; drawpd.h
    ./patchlang &lt; gstate.hh &gt; gstate.h
    ./patchlang &lt; lstsqr.hh &gt; lstsqr.h
    ./patchlang &lt; multipoly.hh &gt; multipoly.h
    ./patchlang &lt; stringo.hh &gt; stringo.h
    ./patchlang &lt; array.hh &gt; array.h
    ./patchlang &lt; equil.hh &gt; equil.h
    ./patchlang &lt; integer.hh &gt; integer.h
    ./patchlang &lt; machdep.hh &gt; machdep.h
    ./patchlang &lt; parse.hh &gt; parse.h
    ./patchlang &lt; teci.hh &gt; teci.h
    ./patchlang &lt; arraylist.hh &gt; arraylist.h
    ./patchlang &lt; fftn.hh &gt; fftn.h
    ./patchlang &lt; keci.hh &gt; keci.h
    ./patchlang &lt; kmeci.hh &gt; kmeci.h
    ./patchlang &lt; mclib.hh &gt; mclib.h
    ./patchlang &lt; phonlib.hh &gt; phonlib.h
    ./patchlang &lt; vectmac.hh &gt; vectmac.h
    ./patchlang &lt; calccorr.hh &gt; calccorr.h
    ./patchlang &lt; findsym.hh &gt; findsym.h
    ./patchlang &lt; lattype.hh &gt; lattype.h
    ./patchlang &lt; misc.hh &gt; misc.h
    ./patchlang &lt; plugin.hh &gt; plugin.h
    ./patchlang &lt; version.hh &gt; version.h
    ./patchlang &lt; calcmf.hh &gt; calcmf.h
    ./patchlang &lt; fixagg.hh &gt; fixagg.h
    ./patchlang &lt; linalg.hh &gt; linalg.h
    ./patchlang &lt; mmclib.hh &gt; mmclib.h
    ./patchlang &lt; predrs.hh &gt; predrs.h
    ./patchlang &lt; xtalutil.hh &gt; xtalutil.h
    ./patchlang &lt; chull.hh &gt; chull.h
    ./patchlang &lt; fxvector.hh &gt; fxvector.h
    ./patchlang &lt; linklist.hh &gt; linklist.h
    ./patchlang &lt; mrefine.hh &gt; mrefine.h
    ./patchlang &lt; refine.hh &gt; refine.h
    ./patchlang &lt; clus_str.hh &gt; clus_str.h
    ./patchlang &lt; getvalue.hh &gt; getvalue.h
    ./patchlang &lt; linsolve.hh &gt; linsolve.h
    ./patchlang &lt; mteci.hh &gt; mteci.h
    ./patchlang &lt; ridge.hh &gt; ridge.h
    ./patchlang &lt; kspacecs.c++ &gt; kspacecs.cc
    ./patchlang &lt; tensor.hh &gt; tensor.h
    ./patchlang &lt; gceutil.hh &gt; gceutil.h
    ./patchlang &lt; tensorsym.hh &gt; tensorsym.h
    ./patchlang &lt; binstream.hh &gt; binstream.h
    ./patchlang &lt; mpiinterf.hh &gt; mpiinterf.h
    ./patchlang &lt; normal.hh &gt; normal.h
    ./patchlang &lt; apb.hh &gt; apb.h
    ./patchlang &lt; opti.hh &gt; opti.h
    ./patchlang &lt; strinterf.hh &gt; strinterf.h
    ./patchlang &lt; meshutil.hh &gt; meshutil.h
    ./patchlang &lt; headers.c++ &gt; headers.cc
    g++ -ffriend-injection -O3    -c -o headers.o headers.cc
    ./patchlang &lt; stringo.c++ &gt; stringo.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported
    g++: warning: switch ‘-ffriend-injection’ is no longer supported


    g++ -ffriend-injection -O3    -c -o stringo.o stringo.cc
    ./patchlang &lt; parse.c++ &gt; parse.cc
    g++ -ffriend-injection -O3    -c -o parse.o parse.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported
    In file included from parse.h:9,
                     from parse.cc:2:
    clus_str.h: In member function ‘int StructureBank&lt;T&gt;::find_new_structures(int)’:
    clus_str.h:202:3: warning: no return statement in function returning non-void [-Wreturn-type]
      202 |   }
          |   ^


    ./patchlang &lt; xtalutil.c++ &gt; xtalutil.cc
    g++ -ffriend-injection -O3    -c -o xtalutil.o xtalutil.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported


    ./patchlang &lt; integer.c++ &gt; integer.cc
    g++ -ffriend-injection -O3    -c -o integer.o integer.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported


    ./patchlang &lt; findsym.c++ &gt; findsym.cc
    g++ -ffriend-injection -O3    -c -o findsym.o findsym.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported


    ./patchlang &lt; clus_str.c++ &gt; clus_str.cc
    g++ -ffriend-injection -O3    -c -o clus_str.o clus_str.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported
    In file included from clus_str.cc:1:
    clus_str.h: In member function ‘int StructureBank&lt;T&gt;::find_new_structures(int)’:
    clus_str.h:202:3: warning: no return statement in function returning non-void [-Wreturn-type]
      202 |   }
          |   ^


    ./patchlang &lt; getvalue.c++ &gt; getvalue.cc
    g++ -ffriend-injection -O3    -c -o getvalue.o getvalue.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported
    getvalue.cc: In function ‘int get_atat_root(AutoString*)’:
    getvalue.cc:257:43: warning: control reaches end of non-void function [-Wreturn-type]
      257 |   AutoString configfilename(getenv("HOME"));
          |                                           ^


    ./patchlang &lt; lattype.c++ &gt; lattype.cc
    g++ -ffriend-injection -O3    -c -o lattype.o lattype.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported


    ./patchlang &lt; mpiinterf.c++ &gt; mpiinterf.cc
    g++ -ffriend-injection -O3    -c -o mpiinterf.o mpiinterf.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported


    ./patchlang &lt; genstr.c++ &gt; genstr.cc
    g++ -ffriend-injection -O3    -c -o genstr.o genstr.cc


    g++: warning: switch ‘-ffriend-injection’ is no longer supported
    In file included from genstr.cc:2:
    clus_str.h: In member function ‘int StructureBank&lt;T&gt;::find_new_structures(int)’:
    clus_str.h:202:3: warning: no return statement in function returning non-void [-Wreturn-type]
      202 |   }
          |   ^


    g++   genstr.o headers.o stringo.o parse.o xtalutil.o integer.o findsym.o clus_str.o getvalue.o lattype.o mpiinterf.o  -lm -o genstr
    rm genstr.o
    make[1]: Leaving directory '/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/packages/atat/src'


### Copy genstr to the current directory.


```python
shutil.copy(f'{atatroot}/src/genstr', f'{rootdir}/')
```




    '/home/jovyan/repos/eneos/matlantis-contrib-binary_alloy/matlantis_contrib_examples/Convex_hull_AlCu_binary_alloy_example/genstr'


