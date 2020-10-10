from setuptools import setup

setup(
   name='classomfier',
   version='0.2',
   author='Javier F. Troncoso',
   author_email='javierfdeztroncoso@gmail.com',
   packages=['classomfier', 'classomfier.test','classomfier.bin'],
   package_data={'classomfier.bin':['classomfier.f90','nomatrix.f90']},
   include_package_data=True,
   url='https://github.com/JaviFdezT/ClasSOMfier',
   license='LICENSE.txt',
   description='ClasSOMfier: A neural network for cluster analysis and detection of lattice defects',
   long_description=open('README.md').read(),
   long_description_content_type="text/markdown",
   install_requires=[],
   download_url="https://github.com/JaviFdezT/ClasSOMfier/archive/0.2.tar.gz",
   keywords = ['kohonen', 'neural', 'network', 'cluster', 'analysis'],
)
