#This is a Python-Package for VirtualFlow#

#You can use the package as below.#

#####step:

    1: cd pyvs
    2: python setup.py install
       now the package support python2,python3
    3: change the all.ctrl,todo.all,other.conf in the VFVS_GK
        you can find the rules in https://docs.virtual-flow.org/tutorials/-LdE94b2AVfBFT72zK-v/vfvs-tutorial-1/setting-up-the-workflow
    4: now you can run follow the command
        start with IPython:
        '''
        from pyvs import pyvs_run
        a=pyvs_run.PyVS('path/to/VFVS_GK')
        a.run()
        '''
#####Tips:
    1:The package is now building.
    2:If you meet any problem when you use the package,please contact with me 'luskyqi@outlook.com'
