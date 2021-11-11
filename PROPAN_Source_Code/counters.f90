!-----------------------------------------------------------------------------------------------!
!    Set Counters                                                                               !
!    Copyright (C) 2021  J. Baltazar                                                            !
!                                                                                               !
!    This program is free software: you can redistribute it and/or modify it under the terms of !
!    the GNU Affero General Public License as published by the Free Software Foundation, either !
!    version 3 of the License, or (at your option) any later version.                           !
!                                                                                               !
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;  !
!    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  !
!    See the GNU Affero General Public License for more details.                                !
!                                                                                               !
!    You should have received a copy of the GNU Affero General Public License                   !
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.                     !
!-----------------------------------------------------------------------------------------------!
SUBROUTINE COUNTERS
!-----------------------------------------------------------------------------------------------!
!    Created by: 29102013, J. Baltazar, version 1.0                                             !
!    Modified  : 16112013, J. Baltazar, version 1.0                                             !
!    Modified  : 26052014, J. Baltazar, Wake Alignment Module                                   !
!    Modified  : 09102014, J. Baltazar, version 2.1                                             !
!    Modified  : 03112014, J. Baltazar, version 2.2                                             !
!    Modified  : 10112014, J. Baltazar, version 3.1, Unsteady Cavitation Model                  !
!    Modified  : 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!    Modified  : 10122014, J. Baltazar, version 3.4, Unsteady Super-Cavitation Model            !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I
!-----------------------------------------------------------------------------------------------!
!    Blade                                                                                      !
!-----------------------------------------------------------------------------------------------!
NRP   =IP*(NRP1-1)
NCP   =IP*(NCP1-1)
NC    =IP*(NCP/2)
NC1   =IP*(NC+1)
NPPAN =IP*(NRP*NCP)
!-----------------------------------------------------------------------------------------------!
!    Write Grid Dimensions                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   WRITE(20,*) 'Spanwise number of panels on the blade NRP=',NRP
   WRITE(20,*) 'Chordwise number of panels on the blade NCP=',NCP
   WRITE(20,*) 'Number of panels on the blade NPPAN=',NPPAN
   WRITE(20,*)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Blade Wake                                                                                 !
!-----------------------------------------------------------------------------------------------!
NRW   =IP*(NRW1-1)
NPW   =IP*(NPW1-1)
NPWPAN=IP*NRW*NPW
NPWCAV=IP*NRW*IABS(NCPW)
!-----------------------------------------------------------------------------------------------!
!    Write Grid Dimensions                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IP == 1) THEN
   WRITE(20,*) 'Spanwise number of panels on the blade wake NRW=',NRW
   WRITE(20,*) 'Streamwise number of panels on the blade wake NPW=',NPW
   WRITE(20,*) 'Number of panels on the blade NPWPAN=',NPWPAN
   WRITE(20,*)
END IF !(IP == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle                                                                                     !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   I=1
   DO WHILE (DABS(XN(I,NNTT/2)-XP(NC1,NRP1)) > TOL)
      I=I+1
   END DO
   NNU=I-1
END IF !(IN == 1)
NNT   =IN*(NNTT/2-1)
NNT1  =IN*(NNT+1)
NNT2  =IN*(NNT+2)
NNTP  =IN*(2*NNT)
NNXT1 =IN*(NNXT-1)
NNX   =IN*(NNXT1/2)
NND   =IN*(NNX-NNU-NC)
NNX1  =IN*(NNX+1)
NNTT1 =IN*(NNTT-1)
NNU1  =IN*(NNU+1)
NN2   =IN*(NNU1+NC)
NCN   =IN*NC
NNPAN =IN*(NNXT1*NNTP)
!-----------------------------------------------------------------------------------------------!
!    Write Grid Dimensions                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   WRITE(20,*) 'Circumferential number of panels on nozzle sector NNTP=',NNTP
   WRITE(20,*) 'Chordwise number of panels on the nozzle NNXT1=',NNXT1
   WRITE(20,*) 'Number of panels on nozzle sector NNPAN=',NNPAN
   WRITE(20,*)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Nozzle Wake                                                                                !
!-----------------------------------------------------------------------------------------------!
NNW   =IN*(NNW1-1)
NNWPAN=IN*NNW*NNTP
!-----------------------------------------------------------------------------------------------!
!    Write Grid Dimensions                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IN == 1) THEN
   WRITE(20,*) 'Circumferential number of panels on nozzle wake sector NNTP=',NNTP
   WRITE(20,*) 'Streamwise number of panels on the nozzle wake NNW=',NNW
   WRITE(20,*) 'Number of panels on nozzle wake sector NNWPAN=',NNWPAN
   WRITE(20,*)
END IF !(IN == 1)
!-----------------------------------------------------------------------------------------------!
!    Hub                                                                                        !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   I=1
   DO WHILE (DABS(XH(I,NHTT/2)-XP(NC1,1)) > TOL)
      I=I+1
   END DO
   NHU=I-1
END IF !(IABS(IH) == 1)
NHX   =IABS(IH)*(NHX1-1)
NHT   =IABS(IH)*(NHTT/2-1)
NHT1  =IABS(IH)*(NHT+1)
NHT2  =IABS(IH)*(NHT+2)
NHTP  =IABS(IH)*(2*NHT)
NHTP1 =IABS(IH)*(NHTP+1)
NHTT1 =IABS(IH)*(NHTT-1)
NHU1  =IABS(IH)*(NHU+1)
NH2   =IABS(IH)*(NHU+NC+1)
NHD   =IABS(IH)*(NHX-NHU-NC)
NCH   =IABS(IH)*NC
NHPAN =IABS(IH)*(NHX*NHTP)
!-----------------------------------------------------------------------------------------------!
!    Write Grid Dimensions                                                                      !
!-----------------------------------------------------------------------------------------------!
IF (IABS(IH) == 1) THEN
   WRITE(20,*) 'Circumferential number of panels on hub sector NHTP=',NHTP
   WRITE(20,*) 'Axial number of panels on the hub NHX=',NHX
   WRITE(20,*) 'Number of panels on hub sector NHPAN=',NHPAN
   WRITE(20,*)
END IF !(IABS(IH) == 1)
!-----------------------------------------------------------------------------------------------!
!    Global                                                                                     !
!-----------------------------------------------------------------------------------------------!
NPAN  =NPPAN+NNPAN+NHPAN+NPWCAV
NWPAN =NPWPAN+NNWPAN
NRT   =NRW+NNTP
NPAN1 =NPAN+NRT
!-----------------------------------------------------------------------------------------------!
!    Write Global Grid Dimensions                                                               !
!-----------------------------------------------------------------------------------------------!
WRITE(20,*) 'Total number of body panels NPAN=',NPAN
WRITE(20,*) 'Total number of wake panels NWPAN=',NWPAN
!-----------------------------------------------------------------------------------------------!
!    Cavitation                                                                                 !
!-----------------------------------------------------------------------------------------------!
NPCAV =0
NSCAV =0
NWCAV =0
!-----------------------------------------------------------------------------------------------!
!    Unsteady                                                                                   !
!-----------------------------------------------------------------------------------------------!
NT    =NREV*NTETA
P     =DFLOAT(NTETA)/DFLOAT(NB)
DTETA =0.D0
IF (NTETA /= 0) DTETA=2.D0*PI/DFLOAT(NTETA)
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE COUNTERS
!-----------------------------------------------------------------------------------------------!
