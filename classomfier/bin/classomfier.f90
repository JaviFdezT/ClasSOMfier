
PROGRAM ClasSOMfier

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



  !Definition of variables and parameters.
  use omp_lib
  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit
  IMPLICIT NONE

  CHARACTER(100) :: arg,inputfile,latparchar,outputfile,resultsfile,folder
  CHARACTER(100) :: outputschar,epochschar,lrchar,sigmachar,usemasschar,pbcchar,useenergychar
  INTEGER  :: rows,outputs,epochs,features=38
  REAL  :: latpar,Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax,learning_rate=0.5,sigma=1.0

  REAL , DIMENSION(6) :: cutoffs
  REAL , DIMENSION(4) :: eta
  REAL , DIMENSION(:), allocatable :: atomsid,masses,energy
  REAL , DIMENSION(:,:), allocatable :: positions,distances 
  REAL , DIMENSION(:,:), allocatable :: inputdata

  LOGICAL :: usemass,useexisting,pbc,useenergy,excludeboundary=.false.

  INTEGER :: nb_ticks_initial, nb_ticks_final, nb_ticks_max, nb_ticks_sec, nb_ticks
  REAL :: elapsed_time

  CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
  CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)

  !First, make sure the right number of inputs have been provided. Then, inputs are read.
  IF(COMMAND_ARGUMENT_COUNT().NE.13) THEN
    WRITE(*,*)'ERROR: 13 COMMAND-LINE ARGUMENTS REQUIRED. STOPPING'
    STOP
  ENDIF

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
  else if (pbcchar.eq."2") then 
          pbc=.false.
          excludeboundary=.true.
  else
    WRITE(*,*)'ERROR. INVALID PARAMETER: pbc, STOPPING'
    STOP
  end if
  

  !Implements a Kohonen network as defined through the input parameters.
  if (.not.useexisting) then 
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
    WRITE ( unit=output_unit, fmt='(t20,a)'   ) "Reading Data..."
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
          CALL readdata
          CALL getdistances
    WRITE ( unit=output_unit, fmt='(t17,a)'   ) "Generating Train Set..."
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
          CALL getinputset
          CALL writedata
  else
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
    WRITE ( unit=output_unit, fmt='(t20,a)'   ) "Reading Data..."
    WRITE ( unit=output_unit, fmt='(a)'   )   '######################################################'
          CALL readdata
          CALL genereteddata
  end if
    
  WRITE ( unit=output_unit, fmt='(t19,a)'   ) "Training Model..."
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
    if (pbc) then
     WRITE ( unit=output_unit, fmt='(a,t51,a)'   ) 'Periodic Boundary Conitions:',       "True" 
    else
     WRITE ( unit=output_unit, fmt='(a,t50,a)'   ) 'Periodic Boundary Conitions:',       "False" 
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
        read(0,*)

        allocate(atomsid(rows))
        allocate(masses(rows))
        allocate(positions(rows,3))
        allocate(distances(rows,rows))
        allocate(inputdata(rows,features))
        if (useenergy) then
          allocate(energy(rows))
          allocate(tmparray(6))
        else
          allocate(tmparray(5))
        end if  

        distances=0.0

        DO i=1,rows
           read(0,*) tmparray
           atomsid(i)=tmparray(1)
           masses(i)=tmparray(2)
           positions(i,:)=tmparray(3:)
           if (useenergy) then
              energy(i)=tmparray(6)
           end if  
        END DO

        CLOSE(0)

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

END SUBROUTINE readdata

SUBROUTINE getdistances
        !
        ! !DESCRIPTION:
        ! Calculates distances for all atoms with their neighbours.
        !

        integer  :: i,j
        real  :: distance,dx,dy,dz,dr(3)
        if (pbc) then
         DO i=1,rows
          DO j=1,rows
           if (i>j) then
            dx=positions(i,1)-positions(j,1)
            dr(1)=dx-nint(dx / (Lxmax-Lxmin)) * (Lxmax-Lxmin)
            dy=positions(i,2)-positions(j,2)
            dr(2)=dy-nint(dy / (Lymax-Lymin)) * (Lymax-Lymin)
            dz=positions(i,3)-positions(j,3)
            dr(3)=dz-nint(dz / (Lzmax-Lzmin)) * (Lzmax-Lzmin)
            distance=norm2(dr)
            distances(i,j)=distance
            distances(j,i)=distance
           end if
          END DO
         END DO
        else
         DO i=1,rows
          DO j=1,rows
           if (i>j) then
            dx=positions(i,1)-positions(j,1)
            dr(1)=dx
            dy=positions(i,2)-positions(j,2)
            dr(2)=dy
            dz=positions(i,3)-positions(j,3)
            dr(3)=dz
            distance=norm2(dr)
            distances(i,j)=distance
            distances(j,i)=distance
           end if
          END DO
         END DO
        end if

END SUBROUTINE getdistances


SUBROUTINE getinputset
        !
        ! !DESCRIPTION:
        ! Applies radial functions to the inputs to describe local environments.
        !
       
        integer  :: i,j,k,feature
        real  :: descriptor1,descriptor2,maxin,minin
        real, allocatable :: distundercut1(:),distundercut2(:)
        REAL, PARAMETER :: PI = 3.1415927

        feature=1
        DO i=1,rows
         inputdata(i,feature)=atomsid(i)
         feature=feature+1
         inputdata(i,feature)=masses(i)
         feature=feature+1
         DO j=1,6
          distundercut1=PACK(distances(i,:), (distances(i,:)<cutoffs(j)).and. &
                           (distances(i,:)>1e-20))
          descriptor1=sum((cos(distundercut1*PI/cutoffs(j))+1)*0.5)
          inputdata(i,feature)=size(distundercut1)
          feature=feature+1
          inputdata(i,feature)=descriptor1
          feature=feature+1
          DO k=1,4
            distundercut2=(cos(distundercut1*PI/cutoffs(j))+1)*0.5*exp(-eta(k)*((distundercut1-0.75*cutoffs(j))**2))
            descriptor2=sum(distundercut2)
            inputdata(i,feature)=descriptor2
            feature=feature+1
          END DO
         END DO
         if (useenergy) then
           inputdata(i,feature)=energy(i)
         end if  
         feature=1
        END DO
        DO i=2,features
          maxin=maxval(inputdata(:,i))
          minin=minval(inputdata(:,i))
          if (maxin.ne.minin) then
            inputdata(:,i)=(inputdata(:,i)-minin)/(maxin-minin)
          end if
        END DO
END SUBROUTINE getinputset

SUBROUTINE writedata
        !
        ! !DESCRIPTION:
        ! Writes the input vectors.
        !
       
        integer  :: i
        real, DIMENSION(5)  :: tmparray

        OPEN (0, file =trim(adjustl(folder)) // "/"//  outputfile)
        
        WRITE ( unit=0, fmt='(i20))') int(rows)

        DO i=1,rows
          WRITE ( unit=0, fmt='(41(f20.10))') inputdata(i,1),inputdata(i,2:)
        END DO
        CLOSE(0)
        
END SUBROUTINE writedata


SUBROUTINE genereteddata
        !
        ! !DESCRIPTION:
        ! Reads input vectors.
        !
       
        integer  :: i
        real, DIMENSION(features)  :: tmparray

        OPEN (0, file = trim(adjustl(folder)) // "/"//outputfile, status='old', action='read')
        
        read(0,*) rows
        

        DO i=1,rows
           read(0,*) inputdata(i,:)
        END DO
        close(0)
END SUBROUTINE genereteddata

SUBROUTINE trainmodel
        !
        ! !DESCRIPTION:
        ! Trains the Kohonen network.
        !
       
        character(len=5) :: file_id
        character(len=500) :: file_name,file_name2
        integer  :: r,i,j,k,l,m,inputs,rowsinside
        real, DIMENSION(:,:), allocatable  :: weights
        integer, DIMENSION(1)  :: minvalue
        real, DIMENSION(rows)  :: minvalueall,categories
        real, DIMENSION(:), allocatable  :: discriminant
        real :: deltaw,category,proximity,distance,meancat

        
        CALL writeinfo

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

        DO r=1,rows
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(r,m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall(r)=discriminant(minvalue(1))
          categories(r)=minvalue(1)
        END DO

        meancat=real(sum(categories))/real(size(categories))
        WRITE ( unit=output_unit, fmt='(t27,a,t42,a)'   ) 'Error','avg. category'
        WRITE ( unit=output_unit, fmt='(a,t25,f10.6,t45,f10.6)'   ) 'Initial:', sum(abs(minvalueall))/size(minvalueall),meancat

       if (excludeboundary) then
        DO i=1,epochs
         rowsinside=0
         DO r=1,rows
          if ((positions(r,1)>Lxmin+1.95*latpar).and.(positions(r,1)<Lxmax-1.95*latpar).and. &
              (positions(r,2)>Lymin+1.95*latpar).and.(positions(r,2)<Lymax-1.95*latpar).and. &
              (positions(r,3)>Lzmin+1.95*latpar).and.(positions(r,3)<Lzmax-1.95*latpar)) then
           DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(r,m:)-weights(j,:))**2))
           END DO
           minvalue=MINLOC(discriminant)
           minvalueall(r)=discriminant(minvalue(1))
           categories(r)=minvalue(1)
           DO j=1,outputs
            DO l=1,inputs
             distance=(real(j-minvalue(1))-nint(real(j-minvalue(1))/real(outputs))*outputs)**2
             proximity=exp(-distance/(2*sigma*decay(real(i))))
             deltaw=learning_rate*decay(real(i))*proximity&
                             *(inputdata(r,m+l-1)-weights(j,l))
             weights(j,l)=weights(j,l)+deltaw
           END DO
          END DO
          rowsinside=rowsinside+1
          end if
         END DO
        END DO
       else
        DO i=1,epochs
         DO r=1,rows
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(r,m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall(r)=discriminant(minvalue(1))
          categories(r)=minvalue(1)
          DO j=1,outputs
            DO l=1,inputs
             distance=(real(j-minvalue(1))-nint(real(j-minvalue(1))/real(outputs))*outputs)**2
             proximity=exp(-distance/(2*sigma*decay(real(i))))
             deltaw=learning_rate*decay(real(i))*proximity&
                             *(inputdata(r,m+l-1)-weights(j,l))
             weights(j,l)=weights(j,l)+deltaw
           END DO
          END DO
         END DO
        END DO
       end if
      

        OPEN (0, file = trim(trim(adjustl(folder)) //  '/' // resultsfile))
       if (excludeboundary) then
        WRITE ( unit=0, fmt='(i0,/)') int(rowsinside)
       else
        WRITE ( unit=0, fmt='(i0,/)') int(rows)
       end if

        DO i=1,outputs
            write(file_id, '(i0)') i
            file_name = trim(adjustl(folder)) //  '/trainset' // trim(adjustl(file_id)) // '.dat' 
            file_name2 = trim(adjustl(folder)) //  '/positions' // trim(adjustl(file_id)) // '.xyz' 
            open(i*100,file = trim(file_name))
            open(i*100+3,file = trim(file_name2))
            WRITE ( unit=i*100+3, fmt='(2(f20.10))') Lxmin,Lxmax
            WRITE ( unit=i*100+3, fmt='(2(f20.10))') Lymin,Lymax
            WRITE ( unit=i*100+3, fmt='(2(f20.10))') Lzmin,Lzmax
        END DO

       if (excludeboundary) then
        DO r=1,rows
          if ((positions(r,1)>Lxmin+1.95*latpar).and.(positions(r,1)<Lxmax-1.95*latpar).and. &
              (positions(r,2)>Lymin+2*latpar).and.(positions(r,2)<Lymax-2*latpar).and. &
              (positions(r,3)>Lzmin+2*latpar).and.(positions(r,3)<Lzmax-2*latpar)) then
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(r,m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall(r)=discriminant(minvalue(1))
          category=minvalue(1)
          WRITE ( unit=0, fmt='(f20.1,3(f20.10))') category,positions(r,:)
          WRITE ( unit=int(category)*100, fmt='(40(f20.10))') inputdata(r,1),inputdata(r,2:)
          if (useenergy) then
           WRITE ( unit=int(category)*100+3, fmt='(40(f20.10))') inputdata(r,1),masses(r),positions(r,:),energy(r)
          else
           WRITE ( unit=int(category)*100+3, fmt='(40(f20.10))') inputdata(r,1),masses(r),positions(r,:)
          end if
         end if
        END DO
      else
        DO r=1,rows
          DO j=1,outputs
            discriminant(j)=sqrt(sum((inputdata(r,m:)-weights(j,:))**2))
          END DO
          minvalue=MINLOC(discriminant)
          minvalueall(r)=discriminant(minvalue(1))
          category=minvalue(1)
          WRITE ( unit=0, fmt='(f20.1,3(f20.10))') category,positions(r,:)
          WRITE ( unit=int(category)*100, fmt='(40(f20.10))') inputdata(r,1),inputdata(r,2:)
          if (useenergy) then
           WRITE ( unit=int(category)*100+3, fmt='(40(f20.10))') inputdata(r,1),masses(r),positions(r,:),energy(r)
          else
           WRITE ( unit=int(category)*100+3, fmt='(40(f20.10))') inputdata(r,1),masses(r),positions(r,:)
          end if
        END DO
       end if
        CLOSE(0)
        DO i=1,outputs
         CLOSE(i*100)
         CLOSE(i*100+3)
        END DO



        meancat=real(sum(categories))/real(size(categories))
        WRITE ( unit=output_unit, fmt='(a,t25,f10.6,t45,f10.6)'   ) 'Final:', sum(abs(minvalueall))/size(minvalueall),meancat
        if (sum(abs(weights))/size(weights)>100) then
          print*,"WARNING: large weights."
        end if

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

END PROGRAM ClasSOMfier


