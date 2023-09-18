import os
import shutil

is_local_int = raw_input("enter '1' for local interaction ")
is_local_int = is_local_int=='1'

if is_local_int:
    op_type = 'real(8)'
    print 'compiling using local interaction'
else:
    op_type = 'integer(1)'
    print 'compiling using transverse-field Ising interaction'

is_complex = raw_input("enter '1' for complex numbers ")
is_complex = is_complex=='1'
if is_complex:
    G_type = 'type(cmatrix)'
    thop_type = 'complex(8)'
    print 'compiling with complex numbers'
else:
    G_type = 'type(rmatrix)'
    thop_type = 'real(8)'
    print 'compiling with real numbers'

filelist = ['dqmc.f90', 'greens_function.f90', 'mc_config.f90', 'module_global.f90', 'measure.f90', 'lapack_interface.f90']

for f in filelist:
    shutil.copy(f, f + '.bak')
    with open(f, 'r') as f2:
        contents = f2.read()

    contents = contents.replace('$op_type$', op_type)
    contents = contents.replace('$G_type$', G_type)
    contents = contents.replace('$thop_type$', thop_type)

    if f=='module_global.f90':
        if is_local_int:
            contents = contents.replace('$is_local_int = 0$', 'is_local_int = 1')
        else:
            contents = contents.replace('$is_local_int = 0$', 'is_local_int = 0')
        if is_complex:
            contents = contents.replace('$is_complex = 1$', 'is_complex = 1')
        else:
            contents = contents.replace('$is_complex = 1$', 'is_complex = 0')

    with open(f, 'w') as f2:
        f2.write(contents)

os.system('make')

for f in filelist:
    shutil.move(f, f +'.made')
    shutil.move(f + '.bak', f)
