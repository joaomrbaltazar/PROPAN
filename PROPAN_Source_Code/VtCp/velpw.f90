!-----------------------------------------------------------------------------------------------!
!    Calculation of the tangential perturbation velocities on the blade wake                    !
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
SUBROUTINE VELPW(TT)
!-----------------------------------------------------------------------------------------------!
!    Created by: 27112014, J. Baltazar, version 3.3, Super-Cavitation Model                     !
!-----------------------------------------------------------------------------------------------!
USE PROPAN_MOD
IMPLICIT NONE
INTEGER :: I,J,TT
!-----------------------------------------------------------------------------------------------!
DO J=1,NRW
   DO I=1,IABS(NCPW)
!-----------------------------------------------------------------------------------------------!
!    u1 Derivative                                                                              !
!-----------------------------------------------------------------------------------------------!
      IF (NCPW < 0) THEN
         IF (I == 1) THEN
            VT1PW(I,J)=A1PW(I,J)*POTPWP(1           ,J,TT)+ &
                       B1PW(I,J)*POTPWP(2           ,J,TT)+ &
                       C1PW(I,J)*POTPWP(3           ,J,TT)
         END IF !(I == 1)
         IF (I == IABS(NCPW)) THEN
            VT1PW(I,J)=A1PW(I,J)*POTPWP(IABS(NCPW)-2,J,TT)+ &
                       B1PW(I,J)*POTPWP(IABS(NCPW)-1,J,TT)+ &
                       C1PW(I,J)*POTPWP(IABS(NCPW)  ,J,TT)
         END IF !(I == IABS(NCPW))
         IF ((I /= 1).AND.(I /= IABS(NCPW))) THEN
            VT1PW(I,J)=A1PW(I,J)*POTPWP(I-1         ,J,TT)+ &
                       B1PW(I,J)*POTPWP(I           ,J,TT)+ &
                       C1PW(I,J)*POTPWP(I+1         ,J,TT)
         END IF !((I /= 1).AND.(I /= IABS(NCPW)))
      ELSEIF (NCPW > 0) THEN
         IF (I == 1) THEN
            VT1PW(I,J)=A1PW(I,J)*POTPWS(1           ,J,TT)+ &
                       B1PW(I,J)*POTPWS(2           ,J,TT)+ &
                       C1PW(I,J)*POTPWS(3           ,J,TT)
         END IF !(I == 1)
         IF (I == IABS(NCPW)) THEN
            VT1PW(I,J)=A1PW(I,J)*POTPWS(IABS(NCPW)-2,J,TT)+ &
                       B1PW(I,J)*POTPWS(IABS(NCPW)-1,J,TT)+ &
                       C1PW(I,J)*POTPWS(IABS(NCPW)  ,J,TT)
         END IF !(I == IABS(NCPW))
         IF ((I /= 1).AND.(I /= IABS(NCPW))) THEN
            VT1PW(I,J)=A1PW(I,J)*POTPWS(I-1         ,J,TT)+ &
                       B1PW(I,J)*POTPWS(I           ,J,TT)+ &
                       C1PW(I,J)*POTPWS(I+1         ,J,TT)
         END IF !((I /= 1).AND.(I /= IABS(NCPW)))
      END IF !(NCPW)
!-----------------------------------------------------------------------------------------------!
   END DO !I=1,IABS(NCPW)
END DO !J=1,NRW
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE VELPW
!-----------------------------------------------------------------------------------------------!
