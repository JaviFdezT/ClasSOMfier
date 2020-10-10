"""
@author: Javier F. Troncoso, October 2020.
        Contact: javierfdeztroncoso@gmail.com
"""


import subprocess
import os
import sys
import glob
import shlex
import pkg_resources
import shutil

class ClasSOMfier:
    """
    ClasSOMfier: A neural network for cluster analysis and detection of lattice defects

    Class that classifies atoms according to their environment.
    Unsupervised training using a 1-dimensional Self Organizing Map (SOM) in Fortran.
    
    Created by Javier F. Troncoso, October 2020.
        Contact: javierfdeztroncoso@gmail.com

    
    USE:
        
        The network and its parameters can be initialized using the following commans:
            >>nn=ClasSOMfier(6.43718,2,"dump1000.file")
        Only 3 parameters are necessary: characteristic length, number of clusters and input file.
        The format of the input file is that provided by the dump command in LAMMPS:
            #compute         peratom all pe/atom
            #dump            dumpid2 all custom 1000 dump*.file id mass x y z c_peratom
        The first command calculates and stores the potential energy per atom.
            
        The network is trained using the following command:
            >>nn.execute()
        The final condigurations are written in ./data (default value) and can be easily read by Ovito.
            
        The final configuration can be postprocessed so that it can be used again to find subcategories
        inside a specific category:
            >>nn.postprocess_output()
            
    """
    
    def __init__(self,latpar,outputs,lammpsoutput,border=False,useexisting=False,epochs=100,
                 learning_rate=0.5,sigma=1.0,traininput="inputdata.dat",
                 trainedoutput="output.xyz",directory="",compilegf=True,pbc=True,
                 usemass=False,useenergy=True,usenomatrix=False):
        """
        Initialize model parameters

        Parameters
        ----------
        latpar : float
            Characteristic distance. Comparable with the lattice parameter.
        outputs : integer
            Number of categories/classes.
        lammpsoutput : string
            File containing atom IDs, masses and positions. 
        useexisting : boolean, optional
            Is True, the train set is not generated and the existing values are used. 
            The default is True.
        epochs : integer, optional
            Number of iterations in the training. The default is 1000.
        learning_rate : float, optional
            Step size in the training. The default is 0.5.
        sigma : float, optional
            Characteristic distance between categories. The default is 1.0.
        traininput : string, optional
            File contaning the train data. The default is "inputdata.dat".
        trainedoutput : string, optional
            File containing the final results: category and positions. The default is "output.xyz".
        directory : string, optional
            Directory in which final results are saved. The default is "", which results in "./data".
        compilegf : boolean, optional
            If True, fortran file is not compiled. The default is True.
        pbc : boolean, optional
            If True, periodic boundary conditions apply. The default is True.
        border : boolean, optional
            If True, atoms at the boundaries are excluded and pbc=False. The default is False.
        usemass : boolean, optional
            If True, masses are considered in the training. The default is False.
        useenergy : boolean, optional
            If True, energies per atom are considered in the training. The default is True.
        usenomatrix : boolean, optional
            If True, no matrices are used to train the neural. As a consecuence, simulations can take a long time.
            The use of large matrices depends on the capacity of the computer. If the number of atoms is too large
              (>11000 approx.), matrices cannot be used.
            The default is True.

        Returns
        -------
        None.

        """
        
        lammpsoutput=os.path.abspath(lammpsoutput)
        if directory=="":
            curdir="/".join(lammpsoutput.split("/")[:-1])
            self.directory=os.path.abspath(curdir+"/data")
        else:    
            self.directory=os.path.abspath(directory)
        self.latpar=latpar
        self.lammpsoutput=lammpsoutput
        self.traininput=traininput
        self.trainedoutput=trainedoutput
        self.epochs=epochs
        self.outputs=outputs
        self.learning_rate=learning_rate
        self.sigma=sigma
        self.compilegf=compilegf
        self.usenomatrix=usenomatrix
    
        self.dash='------------------------------------------------------'
        
        if pbc:
            self.pbc=1
        else:
            self.pbc=0
        if border:
            self.pbc=2
        if usemass:
            self.usemass=1
        else:
            self.usemass=0
        if useenergy:
            self.useenergy=1
        else:
            self.useenergy=0
        if useexisting==True:
            self.useexisting=1
            if not os.path.exists(self.directory+"/"+traininput):
                print(self.dash)
                print("WARNING: Train set has not been generated yet.")
                print(self.dash)
                self.useexisting=0
        else:
            self.useexisting=0
            
    
    def create_directory(self):
        """
        Creates the directory in which final results will be saved.

        Returns
        -------
        None.

        """
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        else:
            print(self.dash)
            print("WARNING: Directory already exists: '{0}'. Some files will be overwritten.".format(self.directory))
            print(self.dash)
        
    def check_compiler(self):
        """
        Checks if the gfortran compiler is installed.

        Returns
        -------
        None.

        """
        if self.compilegf:
           try:
               subprocess.run("gfortran --version", shell=True,check=True,stdout=subprocess.DEVNULL)
           except subprocess.CalledProcessError:
               print(self.dash)
               print("ERROR: gfortran compiler not found")
               print(self.dash)
               sys.exit(0)
            
    def compile_fortran(self):
        """
        Compiles the fortran code.

        Returns
        -------
        None.

        """
        if not self.usenomatrix:
            with open(self.lammpsoutput) as f:
                 atoms = sum(1 for _ in f)-9
            if atoms>10000:
                print(self.dash)
                print("WARNING: The number of atoms is large. ")
                print("WARNING: This may result in a future ERROR after the creation of large matrices.")
                print("WARNING: Consider the use of the option: usenomatrix=True.")
                print(self.dash)
            code=str(pkg_resources.resource_string("classomfier.bin", 'classomfier.f90').decode("utf-8"))
        else:
            code=str(pkg_resources.resource_string("classomfier.bin", 'nomatrix.f90').decode("utf-8"))
        try:
            if self.usenomatrix:
               shutil.copy(self.lammpsoutput,"tmp.dat")
            tmpfile="__tmpfotran.f90"
            filetmp=open(tmpfile,"w")
            filetmp.write(code)
            filetmp.close()
            if os.path.exists(self.directory+"/ClasSOMfier"):
                os.remove(self.directory+"/ClasSOMfier")
            subprocess.run("gfortran -mcmodel=large -Werror=line-truncation -w \
                           -o "+self.directory+"/ClasSOMfier  "+tmpfile, shell=True, check=True)
            if os.path.exists(tmpfile):
                os.remove(tmpfile)
        except subprocess.CalledProcessError:
            print(self.dash)
            print("ERROR: Compilation failure.")
            print(self.dash)
            sys.exit(0)       
    
    def run_fortran(self):
        """
        Runs the fortran code to train the model.

        Returns
        -------
        None.

        """
        if self.usenomatrix:
            print(self.dash)
            print("WARNING: The option 'usenomatrix' is activated. ")
            print("WARNING: Depending on the number of atoms, this run might take a long time.")
            print(self.dash)
        try:
            command=self.directory+"/ClasSOMfier "+str(self.useexisting)+         \
                           " "+str(self.latpar)+" "+str(self.outputs)+" "+       \
                           str(self.epochs)+" "+str(self.learning_rate)+" "+     \
                           str(self.sigma)+" "+self.lammpsoutput+" "+            \
                           self.traininput+" "+self.trainedoutput+" "+           \
                           self.directory+" "+str(self.pbc)+" "+                 \
                           " "+str(self.usemass)+" "+str(self.useenergy)
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(self.dash)
            print("ERROR: Process failed. Please, check if input files exist or their format is right.")
            print(self.dash)
            sys.exit(0)
        if self.usenomatrix:
            if os.path.exists("tmp.dat"):
                    os.remove("tmp.dat")
            if os.path.exists("tmpdes.dat"):
                    os.remove("tmpdes.dat")
               
    def execute(self):
        """
        Creates the train data and trains the model.
        The final results are saved in the trainedoutput file, inside the directory. 
        Additionally, the positions and train data associated with each category are saved
           in separated files starting with positionsN and trainingN, where N is the category.
           

        Returns
        -------
        None.

        """
        self.create_directory()
        self.check_repository()
        self.check_compiler()
        self.compile_fortran()
        self.run_fortran()
        
    def check_repository(self):
        if not os.path.exists(self.lammpsoutput):
               print(self.dash)
               print("ERROR: Input file not found")
               print(self.lammpsoutput,"not found")
               print(self.dash)
               sys.exit(0)
        if not os.path.exists(self.directory):
               print(self.dash)
               print("ERROR: Directory could not be created:")
               print(self.directory)
               print(self.dash)
               sys.exit(0)
        if not os.path.exists(self.directory+"/"+self.traininput) and self.useexisting==1:
               print(self.dash)
               print("ERROR: Train set has not been generated yet.")
               print(self.dash)
               sys.exit(0)
        pass
    
    def postprocess_output(self):
        """
        Reads positionsN and trainingN files and postprocess their format so that they can be used as input
        data in a new run.

        Returns
        -------
        None.

        """
        files=glob.glob(self.directory+"/*")
        trainfiles=[file for file in files if "trainset" in file and "dat" in file]
        positionfiles=[file for file in files if "positions" in file and "xyz" in file]
        if not os.path.exists(self.directory+"/data"):
            os.makedirs(self.directory+"/data")
        for filename in trainfiles:
            rows=sum(1 for line in open(filename))
            file = open(filename, "r")
            contents = file.readlines()
            file.close()
            contents.insert(0, str(rows)+"\n")
            file = open(filename, "w")
            contents = "".join(contents)
            file.write(contents)
            file.close()
            shutil.copy(filename,self.directory+"/data/_"+filename.split("/")[-1])
        for filename in positionfiles:
            rows=sum(1 for line in open(filename))-3
            file = open(filename, "r")
            contents = file.readlines()
            file.close()
            contents.insert(3, "\n")
            contents.insert(0, "\n")
            contents.insert(0, str(rows)+"\n")
            contents.insert(0, "\n")
            contents.insert(0, "\n")
            contents.insert(0, "\n")
            file = open(filename, "w")
            contents = "".join(contents)
            file.write(contents)
            file.close()

