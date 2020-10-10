
PROGRAM ClasSOMfier_nm

  !
  ! ClasSOMfier: A neural network for cluster analysis and detection of lattice defects
  !
  ! !DESCRIPTION:
  !  Class that classifies atoms according to their environment.
  !  Unsupervised training using a 1-dimensional Self Organizing Map (SOM) in Fortran.
  !
  ! History:
  !  Created by Javier F. Troncoso, October 2020
  !      Contact: javierfdeztroncoso@gmail.com
  !
  ! !USES:
  !  A Kohonen network is implemented in Fortran for cluster analysis and detection of lattice defects.
  !  The ClasSOMfier is implemented as a simple program. The different methods can be easily implemented  
  !      using the Python interface.
  !
  ! !REQUIREMENTS:
  !      gfortran. At least, the compiler should be able to read 132 characters by line.
  !
  ! !PROGRAM:
  !  Matrices are not used. Therefore, unlimited date can be analyzed, but at expenses of long simulation times.
  !


  !Definition of variables and parameters.
  use omp_lib
  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit
  IMPLICIT NONE
  CHARACTER(50) :: arg,inputfile,latparchar,outputfile,resultsfile,folder
  CHARACTER(50) :: outputschar,epochschar,lrchar,sigmachar,usemasschar,pbcchar,useenergychar
  integer  :: rows,outputs,epochs,features=38
  real  :: latpar,Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax,learning_rate=0.5,sigma=1.0

  REAL , DIMENSION(6) :: cutoffs
  REAL , DIMENSION(4) :: eta

  logical :: usemass,useexisting,pbc,useenergy

  integer :: nb_ticks_initial, nb_ticks_final, nb_ticks_max, nb_ticks_sec, nb_ticks
  real :: elapsed_time

  !First, make sure the right number of inputs have been provided. Then, inputs are read.
  IF(COMMAND_ARGUMENT_COUNT().NE.13) THEN
    WRITE(*,*)'ERROR: 13 COMMAND-LINE ARGUMENTS REQUIRED. STOPPING'
    STOP
  ENDIF


  CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
  CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)

  CALL GET_COMMAND_ARGUMENT(1,arg)
  CALL GET_COMMAND_ARGUMENT(2,latparchar)
  CALL GET_COMMAND_ARGUMENT(3,outputschar)
  CALL GET_COMMAND_ARGUMENT(4,epochschar)
  CALL GET_COMMAND_ARGUMENT(5,lrchar)
  CALL GET_COMMAND_ARGUMENT(6,sigmachar)
  CALL GET_COMMAND_ARGUMENT(7,inputfile)
  CALL GET_COMMAND_ARGUMENT(8,outputfile)
  CALL GET_COMMAND_ARGUMENT(9,resultsfile)
  CALL GET_COMMAND_ARGUMENT(10,folder)
  CALL GET_COMMAND_ARGUMENT(11,pbcchar)
  CALL GET_COMMAND_ARGUMENT(12,usemasschar)
  CALL GET_COMMAND_ARGUMENT(13,useenergychar)
  READ(latparchar,*)latpar
  READ(outputschar,*)outputs
  READ(epochschar,*)epochs
  READ(lrchar,*)learning_rate
  READ(sigmachar,*)sigma

  if (arg.eq."1") then 
          useexisting=.true.
  else if (arg.eq."0") then 
          useexisting=.false.
  else
    WRITE(*,*)'ERROR. INVALID PARAMETER: useexisting, STOPPING'
    STOP
  end if
  if (usemasschar.eq."1") then 
          usemass=.true.
  else if (usemasschar.eq."0") then 
          usemass=.false.
  else
    WRITE(*,*)'ERROR. INVALID PARAMETER: usemass, STOPPING'
    STOP
  end if
  if (useenergychar.eq."1") then 
          useenergy=.true.
          features=features+1
  else if (useenergychar.eq."0") then 
          useenergy=.false.
  else
    WRITE(*,*)'ERROR. INVALID PARAMETER: useenergy, STOPPING'
    STOP
  end if
  if (pbcchar.eq."1") then 
          pbc=.true.
  else if (pbcchar.eq."0") then 
          pbc=.false.
  else
    WRITE(*,*)'ERROR. INVALID PARAMETER: pbc, STOPPING'
    STOP
  end if
  
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
    WRITE ( unit=output_unit, fmt='(t20,a)'   ) "Reading Data..."
          CALL readdata
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
    WRITE ( unit=output_unit, fmt='(t17,a)'   ) "Generating Train Set..."
          CALL getinputset
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'

    
  WRITE ( unit=output_unit, fmt='(t19,a)'   ) "Training Model..."
  CALL writeinfo 
  CALL trainmodel


    CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
  	nb_ticks = nb_ticks_final - nb_ticks_initial
  	if (nb_ticks_final < nb_ticks_initial) then
          nb_ticks = nb_ticks + nb_ticks_max
	end if
  	elapsed_time   = (REAL(nb_ticks) / nb_ticks_sec)/60
  	
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
    WRITE ( unit=output_unit, fmt='(a,t44,f11.5)'   ) "Time (min)= ",elapsed_time
    WRITE ( unit=output_unit, fmt='(a,/)'   )   '######################################################'


CONTAINS 


SUBROUTINE writeinfo
    !
    ! !DESCRIPTION:
    ! This subroutine prints information about the Kohonen network.
    !
        
    WRITE ( unit=output_unit, fmt='(a)'   )   '------------------------------------------------------'
    WRITE ( unit=output_unit, fmt='(t19,a)'   ) "Model Parameters:"
    WRITE ( unit=output_unit, fmt='(a,t38,i17)'   ) 'Outputs:',          outputs
    WRITE ( unit=output_unit, fmt='(a,t38,f17.5)'   ) 'Lattice Parameter:',          latpar
    WRITE ( unit=output_unit, fmt='(a,t38,i17)'   ) 'Epochs:',          epochs
    WRITE ( unit=output_unit, fmt='(a,t38,i17)'   ) 'Atoms:',          rows
    WRITE ( unit=output_unit, fmt='(a,t38,f17.5)'   ) 'Learning Rate:',          learning_rate
    WRITE ( unit=output_unit, fmt='(a,t38,f17.5)'   ) 'Sigma:',          sigma
    WRITE ( unit=output_unit, fmt='(a,t46,a)'   ) 'Use for large systems:',       "Activated" 
    if (pbc) then
     WRITE ( unit=output_unit, fmt='(a,t51,a)'   ) 'Periodic Boundary Conditions:',       "True" 
    else
     WRITE ( unit=output_unit, fmt='(a,t50,a)'   ) 'Periodic Boundary Conditions:',       "False" 
    end if
    WRITE ( unit=output_unit, fmt='(a,t20,a)'   ) 'Input File:',          inputfile
    WRITE ( unit=output_unit, fmt='(a,t20,a)'   ) 'Train Set:',          trim(adjustl(folder)) // "/"//  outputfile 
    WRITE ( unit=output_unit, fmt='(a,t20,a)'   ) 'Output File:',        trim(adjustl(folder)) // "/"//  resultsfile  

  WRITE ( unit=output_unit, fmt='(a)'   )   '------------------------------------------------------'
END SUBROUTINE writeinfo

SUBROUTINE readdata
        !
        ! !DESCRIPTION:
        ! Reads input data from the input file.
        !
       
        integer  :: i
        real, DIMENSION(:), allocatable  :: tmparray

        OPEN (0, file = trim(adjustl(inputfile)), status='old', action='read')
        
        read(0,*)
        read(0,*)
        read(0,*)
        read(0,*) rows
        read(0,*)
        read(0,*)Lxmin,Lxmax
        read(0,*)Lymin,Lymax
        read(0,*)Lzmin,Lzmax

        cutoffs(1)=latpar*0.55
        cutoffs(2)=latpar*0.75
        cutoffs(3)=latpar*1.10
        cutoffs(4)=latpar*1.40
        cutoffs(5)=latpar*1.60
        cutoffs(6)=latpar*1.90
        eta(1)=1.0
        eta(2)=0.5
        eta(3)=0.1
        eta(4)=0.01

        close(0)

END SUBROUTINE readdata



SUBROUTINE getinputset
        !
        ! !DESCRIPTION:
        ! Applies radial functions to the inputs to describe local environments.
        !
       
        integer  :: i,j,k,l,f,feature
        REAL, PARAMETER :: PI = 3.1415927
        real  :: distance,dx,dy,dz,dr(3),radfunc
        real, DIMENSION(:), allocatable  :: tmparray1,tmparray2
        REAL :: atom,mass,energy
        REAL , DIMENSION(3) :: positions1,positions2
        REAL , DIMENSION(features) :: inputdata,maxin,minin
 
        maxin=-1e-30
        minin= 1e-30
        if (useenergy) then
          allocate(tmparray1(6))
          allocate(tmparray2(6))
        else
          allocate(tmparray1(5))
          allocate(tmparray2(5))
        end if  


        OPEN (1, file = trim(adjustl(inputfile)), status='old', action='read')
        OPEN (200+rows, file ="tmpdes.dat")
        
        WRITE ( unit=200+rows, fmt='(i20))') int(rows)

        
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) 
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)

        OPEN (99, file = "tmp.dat", status='old', action='read')
        DO i=1,rows
          inputdata=0.0
           read(1,*) tmparray1
           atom=tmparray1(1)
           mass=tmparray1(2)
           positions1(:)=tmparray1(3:)
           if (useenergy) then
              energy=tmparray1(6)
           end if  

           inputdata(1)=atom
           inputdata(2)=mass

           read(99,*)
           read(99,*)
           read(99,*)
           read(99,*) 
           read(99,*)
           read(99,*)
           read(99,*)
           read(99,*)
           read(99,*)
         DO j=1,rows

           read(99,*) tmparray2
           positions2(:)=tmparray2(3:)

          if (i.ne.j) then
            dx=positions1(1)-positions2(1)
            dr(1)=dx-nint(dx / (Lxmax-Lxmin)) * (Lxmax-Lxmin)
            dy=positions1(2)-positions2(2)
            dr(2)=dy-nint(dy / (Lymax-Lymin)) * (Lymax-Lymin)
            dz=positions1(3)-positions2(3)
            dr(3)=dz-nint(dz / (Lzmax-Lzmin)) * (Lzmax-Lzmin)
            distance=norm2(dr)
            
            feature=3
            DO k=1,6
              if (distance<cutoffs(k)) then
                radfunc=(cos(distance*PI/cutoffs(k))+1)*0.5
                inputdata(feature)=inputdata(feature)+1
                feature=feature+1
                inputdata(feature)=inputdata(feature)+radfunc
                feature=feature+1
                DO l=1,4
                   inputdata(feature)=inputdata(feature)+radfunc*exp(-eta(l)*((distance-0.75*cutoffs(k))**2))
                   feature=feature+1
                END DO
              end if
            END DO
          end if

         END DO
         if (useenergy) then
           inputdata(features)=energy
         end if  

          rewind(99)
          WRITE ( unit=200+rows, fmt='(41(f20.10))') inputdata(:)
          do f=3,features
             if (inputdata(f)>maxin(f)) then
                 maxin(f)= inputdata(f)
             end if
             if (inputdata(f)<minin(f)) then
                 minin(f)= inputdata(f)
             end if
          end do
        END DO
        close(99)

        close(1)
    rewind(200+rows)

        OPEN (300+rows, file =trim(adjustl(folder)) // "/"//  outputfile)
        WRITE ( unit=300+rows, fmt='(i20))') int(rows)
        read(200+rows,*)
         DO j=1,rows
           read(200+rows,*) inputdata
           do f=3,features
              if (maxin(f).ne.minin(f)) then
                inputdata(f)=(inputdata(f)-minin(f))/(maxin(f)-minin(f))
              end if
           end do
           WRITE ( unit=300+rows, fmt='(41(f20.10))') inputdata(:)
         END DO
        CLOSE(200+rows)
        CLOSE(300+rows)

END SUBROUTINE getinputset



SUBROUTINE trainmodel
        !
        ! !DESCRIPTION:
        ! Trains the Kohonen network.
        !
              
        character(len=5) :: file_id
        character(len=500) :: file_name,file_name2
        integer  :: r,i,j,k,l,m,inputs
        real, DIMENSION(:,:), allocatable  :: weights
        integer, DIMENSION(1)  :: minvalue
        real  :: minvalueall,categories
        real, DIMENSION(:), allocatable  :: discriminant
        real :: deltaw,category,proximity,distance
        REAL , DIMENSION(features) :: inputdata
        real, DIMENSION(:), allocatable  :: tmparray
        REAL , DIMENSION(3) :: positions
        


        if (useenergy) then
          allocate(tmparray(6))
        else
          allocate(tmparray(5))
        end if  
        if (usemass) then
           m=2
           inputs=features-1
        else
           m=3
           inputs=features-2
        end if
        allocate(weights(outputs,inputs))
        allocate(discriminant(outputs))
        call RANDOM_NUMBER(weights)
        weights=(2*weights-1)

        
        OPEN (400+rows, file =trim(adjustl(folder)) // "/"//  outputfile)
        read(400+rows,*)
        minvalueall=0.0
        categories=0.0
        DO r=1,rows
          read(400+rows,*) inputdata
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall=minvalueall+discriminant(minvalue(1))
          categories=categories+minvalue(1)
        END DO

        WRITE ( unit=output_unit, fmt='(t27,a,t42,a)'   ) 'Error','avg. category'
        WRITE ( unit=output_unit, fmt='(a,t25,f10.6,t45,f10.6)'   ) 'Initial:', minvalueall/rows,categories/rows



        DO i=1,epochs
         rewind(400+rows)
         read(400+rows,*)
         DO r=1,rows
          read(400+rows,*) inputdata
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall=minvalueall+discriminant(minvalue(1))
          categories=categories+minvalue(1)
          DO j=1,outputs
            DO l=1,inputs
             distance=(real(j-minvalue(1))-nint(real(j-minvalue(1))/real(outputs))*outputs)**2
             proximity=exp(-distance/(2*sigma*decay(real(i))))
             deltaw=learning_rate*decay(real(i))*proximity&
                             *(inputdata(m+l-1)-weights(j,l))
             weights(j,l)=weights(j,l)+deltaw
           END DO
          END DO
         END DO
        END DO
      


        OPEN (500+rows, file = trim(trim(adjustl(folder)) //  '/' // resultsfile))
        WRITE ( unit=500+rows, fmt='(i0,/)') int(rows)
        OPEN (600+rows, file = trim(adjustl(inputfile)), status='old', action='read')
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)
        read(600+rows,*)


        rewind(400+rows)
        read(400+rows,*)
        minvalueall=0.0
        categories=0.0
        DO r=1,rows
          read(400+rows,*) inputdata
          read(600+rows,*) tmparray
          positions(:)=tmparray(3:)
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall=minvalueall+discriminant(minvalue(1))
          categories=categories+minvalue(1)
          WRITE ( unit=500+rows, fmt='(f20.1,3(f20.10))') real(minvalue(1)),positions(:)
        END DO

        WRITE ( unit=output_unit, fmt='(a,t25,f10.6,t45,f10.6)'   ) 'Final:', minvalueall/rows,categories/rows
        if (sum(abs(weights))/size(weights)>100) then
          print*,"WARNING: large weights."
        end if

        CLOSE(400+rows)
        CLOSE(500+rows)
        CLOSE(600+rows)



END SUBROUTINE trainmodel


function decay(inputval) result(outputval)
    !
    ! !DESCRIPTION:
    ! This function describes the dependence of the learning rate and neighborhood function on the current iteration, inputval.
    !
    real, intent(in) :: inputval ! input
    REAL            :: outputval ! output
    outputval=1./inputval
END FUNCTION decay

END PROGRAM ClasSOMfier_nm


