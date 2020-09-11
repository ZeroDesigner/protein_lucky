**This is a Python-Package for VirtualFlow


*You can use the package as below.*

You can use the package as described below.

Prepare three files:all.ctrl,todo.all,other.conf

Support:python2,python3

**Step:

    1:Prepare three files:all.ctrl,todo.all,other.conf
    2:The other.conf contains the information about './vf_start_jobline.sh 1 12 templates/template1.slurm.sh submit 1 '
    3:Download the package
    4:Install the package with the command lines below
      python setup.py install
    5:Change the all.ctrl,todo.all,other.conf in the VFVS_GK
      you can find the rules in https://docs.virtual-flow.org/tutorials/-LdE94b2AVfBFT72zK-v/vfvs-tutorial-1/setting-up-the-workflow
    4:Now you can run follow the command
        start with IPython:
        '''
        from pyvs import pyvs_run
        a=pyvs_run.PyVS('path/to/VFVS_GK')
        a.run()
        '''
**Tips:
    1:The package is now building.
    2:If you meet any problem when you use the package,please contact with me 'luskyqi@outlook.com'
