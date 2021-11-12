!-----------------------------------------------------------------------------------------------!
!    Fourier Coefficients                                                                       !
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
SUBROUTINE FOURIER_COEF(NI,TI,XT,FT,AN,BN)
!-----------------------------------------------------------------------------------------------!
!    Created by: Joao Baltazar, IST                                                             !
!    Modified  : 22062016, J. Baltazar, non-uniform spacing                                     !
!-----------------------------------------------------------------------------------------------!
IMPLICIT NONE
INTEGER :: I,IM1,J,NI,TI
DOUBLE PRECISION :: PI,XT(TI),FT(TI),MCOS(TI,NI),MSIN(TI,NI)
DOUBLE PRECISION :: AN(0:NI),BN(1:NI)
!-----------------------------------------------------------------------------------------------!
PI=4.D0*DATAN(1.D0)
!-----------------------------------------------------------------------------------------------!
!    Cosine Matrix                                                                              !
!-----------------------------------------------------------------------------------------------!
DO J=1,NI
   DO I=1,TI
      MCOS(I,J)=DCOS(DFLOAT(J)*XT(I))
   END DO !I=1,TI
END DO !J=1,NI
!-----------------------------------------------------------------------------------------------!
!    Sine Matrix                                                                                !
!-----------------------------------------------------------------------------------------------!
DO J=1,NI
   DO I=1,TI
      MSIN(I,J)=DSIN(DFLOAT(J)*XT(I))
   END DO !I=1,TI
END DO !J=1,NI
!-----------------------------------------------------------------------------------------------!
!    A0 Coefficient                                                                             !
!-----------------------------------------------------------------------------------------------!
AN(0)=0.D0
DO I=2,TI
   IM1=I-1
   AN(0)=AN(0)+(FT(I)+FT(IM1))*0.5D0*(XT(I)-XT(IM1))
END DO !I=2,TI
AN(0)=AN(0)/PI
!-----------------------------------------------------------------------------------------------!
!    An Coefficients                                                                            !
!-----------------------------------------------------------------------------------------------!
DO J=1,NI
   AN(J)=0.D0
   DO I=2,TI
      IM1=I-1
      AN(J)=AN(J)+(FT(I)*MCOS(I,J)+FT(IM1)*MCOS(IM1,J))*0.5D0*(XT(I)-XT(IM1))
   END DO !I=2,TI
   AN(J)=AN(J)/PI
   IF (DABS(AN(J)) < 1.D-14) AN(J)=0.D0
END DO !J=1,NI
!-----------------------------------------------------------------------------------------------!
!    Bn Coefficients                                                                            !
!-----------------------------------------------------------------------------------------------!
DO J=1,NI
   BN(J)=0.D0
   DO I=2,TI
      IM1=I-1
      BN(J)=BN(J)+(FT(I)*MSIN(I,J)+FT(IM1)*MSIN(IM1,J))*0.5D0*(XT(I)-XT(IM1))
   END DO !I=2,TI
   BN(J)=BN(J)/PI
   IF (DABS(BN(J)) < 1.D-14) BN(J)=0.D0
END DO !J=1,NI
!-----------------------------------------------------------------------------------------------!
END SUBROUTINE FOURIER_COEF
!-----------------------------------------------------------------------------------------------!
