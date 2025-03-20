PROGRAM traject_vtf

!***********************************************************************
!
!     daresbury laboratory ccp5 program for post-processing HISTORY
!     files created from dl_meso_dpd to produce trajectory file in
!     vtf format
!
!     copyright ukri stfc daresbury laboratory
!     author - m. a. seaton december 2018
!     amended - m. a. seaton december 2019
!
!***********************************************************************

      IMPLICIT none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND (15, 307)
      INTEGER, PARAMETER :: li = SELECTED_INT_KIND (12)
      INTEGER, PARAMETER :: nuser=5,nrtout=9,ntraj=11
      INTEGER, PARAMETER :: endversion = 1

      CHARACTER(80) :: text
      CHARACTER(8), ALLOCATABLE :: namspe (:), nammol (:)

      INTEGER, ALLOCATABLE :: ltp (:), ltm (:), mole (:), bndtbl (:,:), nspe (:), readint (:), globindex (:)
      INTEGER :: chain, ibond, imxspe, ioerror, i, k, nmoldef, cmdlinarg
      INTEGER :: maxnum, numspe, nbeads, nubeads, numbeads, numbonds, global, species, molecule
      INTEGER :: lfrzn, keytrj, srfx, srfy, srfz, endver, Dlen, numframe, nstep, framesize, lend, leni
      INTEGER(KIND=li) :: filesize, mypos, currentpos, headerpos, lend_li, leni_li, numbeadsli, framesizeli

      REAL(KIND=dp), ALLOCATABLE :: amass (:), bbb (:), chge (:), readdata (:)
      REAL(KIND=dp) :: dimx, dimy, dimz, halfx, halfy, halfz, shrdx, shrdy, shrdz
      REAL(KIND=dp) :: ndx, ndy, ndz, time

      LOGICAL :: volchange, eof, separate, struct, swapend, bigend, sorted

!     determine number of bytes for selected double precision kind
!     (the default SELECTED_REAL_KIND (15, 307) should return 8 bytes)

      lend = STORAGE_SIZE (1.0_dp) / 8
      leni = BIT_SIZE (1) / 8
      lend_li = INT (lend, KIND=li)
      leni_li = INT (leni, KIND=li)

!     check endianness of machine

      bigend = (IACHAR(TRANSFER(1,"a"))==0)

!     determine if separate files for unbonded and bonded beads are requested,
!     separate files for structure and coordinates are requested, particles
!     are to be sorted by index (or not), or if help requested by user

      separate = .false.
      struct = .false.
      sorted = .false.
      cmdlinarg = command_argument_count ()
      IF (cmdlinarg>0) THEN
        DO i = 1, cmdlinarg
          CALL get_command_argument (i, text)
          IF (text (1:2) == '-b' .OR. text (1:2) == '-B') separate = .true.
          IF (text (1:3) == '-sc' .OR. text (1:3) == '-SC') struct = .true.
          IF (text (1:2) == '-s' .OR. text (1:2) == '-S') sorted = .true.
          IF (text (1:2) == '-u' .OR. text (1:2) == '-U') sorted = .false.
          IF (text (1:2) == '-h' .OR. text (1:2) == '-H') THEN
            PRINT *, "DL_MESO traject utility, converts trajectory data in HISTORY file to"
            PRINT *, "VTF format for VMD (default: creates traject.vtf file)"
            CALL get_command_argument (0, text)
            PRINT *, "Usage: ", TRIM(text), " [OPTIONS]"
            PRINT *, " "
            PRINT *, "Options:"
            PRINT *, " "
            PRINT *, "-h"
            PRINT *, "       display this help and exit"
            PRINT *, " "
            PRINT *, "-b"
            PRINT *, "       write separate files for bonded and unbonded particles"
            PRINT *, " "
            PRINT *, "-sc"
            PRINT *, "       write separate files for structure (VSF) and coordinates (VCF)"
            PRINT *, " "
            PRINT *, "-s"
            PRINT *, "       sort particles by global index and write to trajectory file"
            PRINT *, "       in that order"
            PRINT *, " "
            PRINT *, "-u"
            PRINT *, "       write particles to trajectory file in the order they appear"
            PRINT *, "       in HISTORY file (default)"
            PRINT *, " "
            STOP
          END IF
        END DO
      END IF

!     determine if HISTORY file exists, which endianness to use,
!     if type of real is correct

      INQUIRE (file = 'HISTORY', EXIST = eof)
      IF (.NOT. eof) THEN
        PRINT *, "ERROR: cannot find HISTORY file"
        STOP
      END IF

      OPEN (ntraj, file = 'HISTORY', access = 'stream', form = 'unformatted', status = 'unknown')

      swapend = .false.
      READ (ntraj) endver, Dlen

      IF (endver/=endversion) THEN
        swapend = .true.
        CLOSE (ntraj)
        IF (bigend) THEN
          OPEN (ntraj, file = 'HISTORY', access = 'stream', form = 'unformatted', status = 'unknown', convert = 'little_endian')
        ELSE
          OPEN (ntraj, file = 'HISTORY', access = 'stream', form = 'unformatted', status = 'unknown', convert = 'big_endian')
        END IF
        READ (ntraj) endver, Dlen
        IF (endver/=endversion) THEN
          PRINT *, "ERROR: corrupted HISTORY file or created with incorrect version of DL_MESO"
          STOP
        END IF
      END IF

      IF (Dlen/=lend) THEN
        PRINT *, "ERROR: incorrect type of real number used in HISTORY file"
        PRINT *, "       recompile traject.f90 with reals of ", Dlen, " bytes"
        STOP
      END IF

!     read file size, number of frames and timestep numbers

      READ (ntraj) filesize, numframe, nstep

!     read HISTORY file to determine species names and properties,
!     find numbers of beads and bonds and allocate arrays ready for
!     bead, bond and available trajectory data

      READ (ntraj) text

      PRINT *, 'writing trajectory data from HISTORY file for simulation:'
      PRINT *, text

      READ (ntraj) numspe, nmoldef, nubeads, numbeads, numbonds, keytrj, srfx, srfy, srfz

      ALLOCATE (namspe (numspe), nammol (nmoldef))
      ALLOCATE (amass (numspe), bbb (numspe), chge (numspe), ltp (1:numbeads), ltm (1:numbeads), mole (1:numbeads), nspe (numspe))

      DO i = 1, numspe
        READ (ntraj) namspe (i), amass (i), bbb (i), chge (i), lfrzn
      END DO

      DO i = 1, nmoldef
        READ (ntraj) nammol (i)
      END DO

      framesize = (keytrj+1) * 3
      ALLOCATE (bndtbl (numbonds, 2), readdata (1:framesize), readint (1:numbeads), globindex (1:numbeads))

!     if no bonds exist, switch off separate .vtf file option

      IF (numbonds==0) separate = .false.

!     read in data for beads and bonds

      ibond = 0
      nspe = 0

      DO i = 1, numbeads
        READ (ntraj) global, species, molecule, chain
        ltp (global) = species
        ltm (global) = molecule
        mole (global) = chain
        nspe (species) = nspe (species) + 1
      END DO

      IF (numbonds>0) THEN
        DO i = 1, numbonds
          READ (ntraj) bndtbl (i, 1), bndtbl (i, 2)
        END DO
      END IF

!     reached end of header: find current position in file

      INQUIRE (unit=ntraj, POS=headerpos)

!     open output file(s) (.VTF or .VSF for VMD) to write
!     bead types, bonds and box size (if using .VSF file)

      IF (separate) THEN
        IF (struct) THEN
          OPEN (nrtout, file = 'traject_bead.vsf')
          OPEN (nrtout+1, file = 'traject_mole.vsf')
        ELSE
          OPEN (nrtout, file = 'traject_bead.vtf')
          OPEN (nrtout+1, file = 'traject_mole.vtf')
        END IF
      ELSE
        IF (struct) THEN
          OPEN (nrtout, file = 'traject.vsf')
        ELSE
          OPEN (nrtout, file = 'traject.vtf')
        END IF
      END IF

      maxnum = 0
      DO i = 1, numspe
        IF (nspe (i)>maxnum) THEN
          maxnum = nspe (i)
          imxspe = i
        END IF
      END DO

!     define beads

      IF (numbeads==maxnum .AND. numbonds==0) THEN
        WRITE (nrtout, '("atom 0:",I10.10,"    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8)') &
              &numbeads-1, bbb (imxspe), amass (imxspe), chge (imxspe), namspe (imxspe)
      ELSE
        WRITE (nrtout, '("atom default    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8)') &
              &bbb (imxspe), amass (imxspe), chge (imxspe), namspe (imxspe)
        IF (separate) THEN
          WRITE (nrtout+1, '("atom default    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8)') &
              &bbb (imxspe), amass (imxspe), chge (imxspe), namspe (imxspe)
          DO i = 1, numbeads
            IF ((ltp (i)/=imxspe .OR. i==nubeads) .AND. ltm(i)==0) THEN
              WRITE (nrtout, '("atom ",I10,"    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8)') &
                    &i-1, bbb (ltp (i)), amass (ltp (i)), chge (ltp (i)), namspe (ltp (i))
            ELSE IF (ltm (i)/=0) THEN
              WRITE (nrtout+1, '("atom ",I10,"    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8,&
                    &" resid ",I10," resname ",A8)') &
                    &i-nubeads-1, bbb (ltp (i)), amass (ltp (i)), chge (ltp (i)), namspe (ltp (i)), mole (i), nammol (ltm (i))
            END IF

          END DO
        ELSE
          DO i = 1, numbeads
            IF ((ltp (i)/=imxspe .OR. i==numbeads) .AND. ltm(i)==0) THEN
              WRITE (nrtout, '("atom ",I10,"    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8)') &
                    &i-1, bbb (ltp (i)), amass (ltp (i)), chge (ltp (i)), namspe (ltp (i))
            ELSE IF (ltm (i)/=0) THEN
              WRITE (nrtout, '("atom ",I10,"    radius ",F10.6," mass ",F10.6," charge ",F10.6," name ",A8,&
                    &" resid ",I10," resname ",A8)') &
                    &i-1, bbb (ltp (i)), amass (ltp (i)), chge (ltp (i)), namspe (ltp (i)), mole (i), nammol (ltm (i))
            END IF
          END DO
        END IF
      END IF

!     define bonds

      IF (separate) THEN
        WRITE (nrtout+1, '()')
        DO i = 1, numbonds
          WRITE (nrtout+1, '("bond ",I10,":",I10)') (bndtbl (i, 1)-nubeads-1), (bndtbl (i, 2)-nubeads-1)
        END DO
      ELSE IF (numbonds>0 .AND. numbeads>nubeads) THEN
        WRITE (nrtout, '()')
        DO i = 1, numbonds
          IF ((bndtbl (i, 1)-1)>=0) WRITE (nrtout, '("bond ",I10,":",I10)') (bndtbl (i, 1)-1), (bndtbl (i, 2)-1)
        END DO
      END IF

      DEALLOCATE (ltp, ltm, mole, bndtbl, nspe, namspe, amass, bbb, chge)

!     obtain positions and velocities for all beads at each time step

      dimx = 0.0_dp
      dimy = 0.0_dp
      dimz = 0.0_dp
      eof = .false.
      framesizeli = INT (framesize, KIND=li)
      numbeadsli = INT (numbeads, KIND=li)

      DO k = 1, numframe

        currentpos = headerpos + INT (k-1, KIND=li) * ((7_li + numbeadsli * framesizeli) * lend_li + (1_li + numbeadsli) * leni_li)
        READ (ntraj, POS=currentpos, IOSTAT=ioerror) time, nbeads, ndx, ndy, ndz, shrdx, shrdy, shrdz

        IF (ioerror/=0) THEN
          eof = .true.
          IF (k==0) THEN
            PRINT *, 'ERROR: cannot find trajectory data in HISTORY file'
            STOP
          END IF
          EXIT
        END IF

        volchange = (ndx/=dimx .OR. ndy/=dimy .OR. ndz/=dimz)
        dimx = ndx
        dimy = ndy
        dimz = ndz
        halfx = 0.5_dp * dimx
        halfy = 0.5_dp * dimy
        halfz = 0.5_dp * dimz

        ! if writing separate files for structure and coordinates,
        ! write box dimensions to structure file(s) and close,
        ! then open coordinate file(s) (.VCF)

        IF (struct .AND. k==1) THEN
          WRITE (nrtout, '(/,"pbc ", 3F12.6, " 90 90 90")') dimx, dimy, dimz
          IF (separate) WRITE (nrtout+1, '(/,"pbc ", 3F12.6, " 90 90 90")') dimx, dimy, dimz
          CLOSE (nrtout)
          IF (separate) CLOSE (nrtout+1)
          IF (separate) THEN
            OPEN (nrtout, file = 'traject_bead.vcf')
            OPEN (nrtout+1, file = 'traject_mole.vcf')
          ELSE
            OPEN (nrtout, file = 'traject.vcf')
          END IF
        END IF

        IF (volchange .OR. k==1) THEN
          WRITE (nrtout, '(/,"timestep indexed",/,"pbc ", 3F12.6, " 90 90 90")') dimx, dimy, dimz
          IF (separate) WRITE (nrtout+1, '(/,"timestep indexed",/,"pbc ", 3F12.6, " 90 90 90")') dimx, dimy, dimz
        ELSE
          WRITE (nrtout, '(/, "timestep indexed")')
          IF (separate) WRITE (nrtout+1, '(/, "timestep indexed")')
        END IF

        IF (sorted) THEN

          READ (ntraj) readint (1:nbeads)
          CALL quicksort_integer_indexed (readint, 1, nbeads, globindex)
          IF (separate) THEN
            DO i = 1, nubeads
              global = globindex (i)
              mypos = currentpos + leni_li * (1_li + numbeadsli) + (7_li + INT (global-1, KIND=li) * framesizeli) * lend_li
              READ (ntraj, POS=mypos) readdata (1:framesize)
              ! KaGaSi edit (snippet by Micheal Seaton)
              SELECT CASE(keytrj)
                CASE (0)
                  WRITE (nrtout, '(I10,X,3(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz
                CASE (1)
                  WRITE (nrtout, '(I10,X,6(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:6)
                CASE (2)
                  WRITE (nrtout, '(I10,X,9(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:9)
              END SELECT
              ! KaGaSi edit ends
            END DO
            DO i = nubeads+1, nbeads
              global = globindex (i)
              mypos = currentpos + leni_li * (1_li + numbeadsli) + (7_li + INT (global-1, KIND=li) * framesizeli) * lend_li
              READ (ntraj, POS=mypos) readdata (1:framesize)
              ! KaGaSi edit (snippet by Micheal Seaton)
              SELECT CASE(keytrj)
                CASE (0)
                  WRITE (nrtout, '(I10,X,3(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz
                CASE (1)
                  WRITE (nrtout, '(I10,X,6(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:6)
                CASE (2)
                  WRITE (nrtout, '(I10,X,9(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:9)
              END SELECT
              ! KaGaSi edit ends
            END DO
          ELSE
            DO i = 1, nbeads
              global = globindex (i)
              mypos = currentpos + leni_li * (1_li + numbeadsli) + (7_li + INT (global-1, KIND=li) * framesizeli) * lend_li
              READ (ntraj, POS=mypos) readdata (1:framesize)
              ! KaGaSi edit (snippet by Micheal Seaton)
              SELECT CASE(keytrj)
                CASE (0)
                  WRITE (nrtout, '(I10,X,3(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz
                CASE (1)
                  WRITE (nrtout, '(I10,X,6(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:6)
                CASE (2)
                  WRITE (nrtout, '(I10,X,9(F12.6,X))') i-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:9)
              END SELECT
              ! KaGaSi edit ends
            END DO
          END IF

        ELSE

          READ (ntraj) globindex (1:nbeads)
          IF (separate) THEN
            DO i = 1, nbeads
              READ (ntraj) readdata (1:framesize)
              global = globindex (i)
              IF (global>nubeads) THEN
                ! KaGaSi edit (snippet by Micheal Seaton)
                SELECT CASE(keytrj)
                  CASE (0)
                    WRITE (nrtout, '(I10,X,3(F12.6,X))') global-nubeads-1, readdata(1)+halfx, readdata(2)+halfy, &
                           &readdata(3)+halfz
                  CASE (1)
                    WRITE (nrtout, '(I10,X,6(F12.6,X))') global-nubeads-1, readdata(1)+halfx, readdata(2)+halfy, &
                           &readdata(3)+halfz, readdata(4:6)
                  CASE (2)
                    WRITE (nrtout, '(I10,X,9(F12.6,X))') global-nubeads-1, readdata(1)+halfx, readdata(2)+halfy, &
                           &readdata(3)+halfz, readdata(4:9)
                END SELECT
                ! KaGaSi edit ends
              ELSE
              ! KaGaSi edit (snippet by Micheal Seaton)
              SELECT CASE(keytrj)
                CASE (0)
                  WRITE (nrtout, '(I10,X,3(F12.6,X))') global-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz
                CASE (1)
                  WRITE (nrtout, '(I10,X,6(F12.6,X))') global-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:6)
                CASE (2)
                  WRITE (nrtout, '(I10,X,9(F12.6,X))') global-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:9)
              END SELECT
              ! KaGaSi edit ends
              END IF
            END DO
          ELSE
            DO i = 1, nbeads
              READ (ntraj) readdata (1:framesize)
              ! KaGaSi edit (snippet by Micheal Seaton)
              SELECT CASE(keytrj)
                CASE (0)
                  WRITE (nrtout, '(I10,X,3(F12.6,X))') globindex(i)-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz
                CASE (1)
                  WRITE (nrtout, '(I10,X,6(F12.6,X))') globindex(i)-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:6)
                CASE (2)
                  WRITE (nrtout, '(I10,X,9(F12.6,X))') globindex(i)-1, readdata(1)+halfx, readdata(2)+halfy, readdata(3)+halfz, &
                         &readdata(4:9)
              END SELECT
              ! KaGaSi edit ends
            END DO
          END IF

        END IF

        IF (eof) EXIT

        PRINT *, "Timestep ",k," : time = ",time

      END DO

      CLOSE (ntraj)
      CLOSE (nrtout)
      IF (separate) CLOSE (nrtout+1)

CONTAINS

      SUBROUTINE quicksort_integer_indexed (list, stride, n, indices)

!**********************************************************************
!
!     sort integers in array into numerical order, recording original
!     positions of values (routine to prepare indices array)
!
!     copyright ukri stfc daresbury laboratory
!     authors - m. a. seaton august 2013
!
!**********************************************************************

      INTEGER, INTENT (INOUT) :: list (:)
      INTEGER, INTENT (IN) :: stride, n
      INTEGER, INTENT (OUT) :: indices (:)
      INTEGER :: i

      DO i = 1, n
        indices (i) = i
      END DO

      CALL qsort_integer (list, indices, stride, 1, n)

      END SUBROUTINE quicksort_integer_indexed

      RECURSIVE SUBROUTINE qsort_integer (list, index, stride, low, high)

!**********************************************************************
!
!     sort integers in array into numerical order, recording original
!     positions of values
!
!     copyright ukri stfc daresbury laboratory
!     authors - m. a. seaton august 2013
!
!**********************************************************************

      INTEGER, INTENT (INOUT) :: list (:), index (:)
      INTEGER, INTENT (IN) :: low, high
      INTEGER, INTENT (IN) :: stride
      INTEGER :: i, j, k, reference, temp

      IF (high < low + 6) THEN

!     resort to bubble sort for very small lists (5 items or fewer)

        DO i = low, high - 1
          DO j = i+1, high
            IF (list (stride * (i - 1) + 1) > list (stride * (j - 1) + 1)) THEN
              DO k = 1, stride
                temp = list (stride * (i - 1) + k)
                list (stride * (i - 1) + k) = list (stride * (j - 1) + k)
                list (stride * (j - 1) + k) = temp
              END DO
              temp = index (i)
              index (i) = index (j)
              index (j) = temp
            END IF
          END DO
        END DO

      ELSE

!     apply partition-based sort

        reference = list (stride * ((low+high)/2 - 1) + 1)
        i = low - 1
        j = high + 1
        DO
          DO
            i = i + 1
            IF (list (stride * (i-1) + 1) >= reference) EXIT
          END DO
          DO
            j = j - 1
            IF (list (stride * (j-1) + 1) <= reference) EXIT
          END DO
          IF (i < j) THEN
            DO k = 1, stride
              temp = list (stride * (i-1) + k)
              list (stride * (i-1) + k) = list (stride * (j-1) + k)
              list (stride * (j-1) + k) = temp
            END DO
            temp = index (i)
            index (i) = index (j)
            index (j) = temp
          ELSE IF (i == j) THEN
            i = i + 1
            EXIT
          ELSE
            EXIT
          END IF
        END DO

        IF (low<j) CALL qsort_integer (list, index, stride, low, j)
        IF (i<high) CALL qsort_integer (list, index, stride, i, high)

      END IF

      END SUBROUTINE qsort_integer

END PROGRAM traject_vtf
