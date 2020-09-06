import os
import sys
import shutil
import stat
import zipfile
class PyVS():
    def __init__(self,conf_folder):
        self.conf_folder=conf_folder
    def run(self):
            global tar_path
            tar_path=self.conf_folder+'/VFVS_tools'
            abs_path=os.path.abspath(__file__)
            dir_path=abs_path[:-12]
            dock_bin=dir_path+'/VFVS_tools/bin'
#            VFTools_bin=dir_path+'/VFTools_master/bin'

            os.system('export PATH='+dock_bin+':$PATH')
#            os.system('export PATH='+VFTools_bin+':$PATH')

            shutil.copy(dir_path+"/VFVS_Tools.zip",self.conf_folder)

            f_zip = zipfile.ZipFile(self.conf_folder+"/VFVS_Tools.zip",'r')
            for file in f_zip.namelist():
                f_zip.extract(file)



            shutil.copy(self.conf_folder+'/all.ctrl',self.conf_folder+'/VFVS_tools/templates')
            shutil.copy(self.conf_folder+'/todo.all',self.conf_folder+'/VFVS_tools/templates')

            os.chdir(self.conf_folder+'/VFVS_tools/')
            os.system('chmod +x *')
            os.system('chmod +x */*')
            os.system('./vf_prepare_folders_new.sh')

            f=open(self.conf_folder+'/other.conf','r')
            lines=f.readlines()
            lines=[i.strip() for i in lines]
            f.close()

            f_sys=open(self.conf_folder+'/all.ctrl','r')
            lines_sys=f_sys.readlines() 
            lines_sys=[i.strip() for i in lines_sys if 'batchsystem' in i]
            system_id=lines_sys[1].split('=')[1].lower()
            os.system('./vf_start_jobline.sh '+ lines[1]+' '+lines[3]+' '+'templates/template1.' + system_id +'.sh submit '+lines[1])
            print('./vf_start_jobline.sh '+ lines[1]+' '+lines[3]+' '+'templates/template1.' + system_id +'.sh submit '+lines[1])
