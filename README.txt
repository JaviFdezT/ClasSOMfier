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
