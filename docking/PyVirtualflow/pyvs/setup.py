
from setuptools import setup, find_packages
# 导入setup函数并传参
setup(name='pyvs',
       version='0.1.0',
       description='Python for VirtualFlow',
       author='sujiaqi',
       author_email='luskyqi@outlook.com',
       packages=find_packages(),
       package_data={'':['*.zip']}
       )

