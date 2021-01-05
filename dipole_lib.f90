!f2py -c -m fortran_dip dipole_lib.f90

SUBROUTINE calc_corr(time_series1, time_series2, N, Nmax, correlation)

		IMPLICIT NONE
		INTEGER, PARAMETER :: dp=KIND(0.0D0)
		REAL(KIND=dp) :: Pi
		REAL(KIND=dp), DIMENSION(3,0:N) :: time_series1, time_series2
		REAL(KIND=dp), DIMENSION(0:Nmax) :: correlation
		INTEGER :: N,I,J,Nmax,NF
		
        !f2py intent(in) :: time_series1, time_series2, Nmax
        !f2py intent(hide), depend(N) :: N = len(time_series1)
        !f2py intent(out) correlation

		correlation=0.0_dp
		IF ( 2*Nmax > N) THEN
			WRITE(6,*) "Error! Dipole data (", N, " frames) too short for", Nmax, " ps correlation function! Exiting!"
			 call EXIT(-1)
		END IF
		WRITE(6,*) "Start to calculate dipdip...."
		N = N
		DO I = 0, Nmax
			NF = N - I
		  DO J = 0, NF
		    correlation(I) = correlation(I) + DOT_PRODUCT(time_series1(:,J),time_series2(:,J+I))
		 ENDDO
			correlation(I) = correlation(I)/REAL(NF + 1,kind=dp)
		ENDDO

END SUBROUTINE

SUBROUTINE calc_two_side_corr(time_series1, time_series2, N, Nmax, correlation)

		IMPLICIT NONE
		INTEGER, PARAMETER :: dp=KIND(0.0D0)
		REAL(KIND=dp), DIMENSION(3,0:N) :: time_series1, time_series2
		REAL(KIND=dp), DIMENSION(0:Nmax) :: correlation
		INTEGER :: N,I,J,Nmax

        !f2py intent(in) :: time_series1, time_series2, Nmax
        !f2py intent(hide), depend(N) :: N = len(time_series1)
        !f2py intent(out) correlation
		IF ( 2*Nmax > N) THEN
			WRITE(6,*) "Error! Dipole data (", N, " frames) too short for", Nmax, " ps correlation function! Exiting!"
			 call EXIT(-1)
		END IF

		correlation=0.0_dp
        WRITE(6,*) "Start to calculate dipdip...."
        DO I=1,N-Nmax
            DO J=I,I+Nmax
                correlation(J-I)=correlation(J-I)+DOT_PRODUCT(time_series1(:,I),time_series2(:,J))
            ENDDO
        ENDDO

        DO I=0,Nmax
            correlation(I)=correlation(I)/(REAL(N-I,kind=dp)*REAL(N,kind=dp))
        ENDDO

END SUBROUTINE

SUBROUTINE fl_transform(correlation, time, freq, n_time, n_freq, f_image)
		IMPLICIT NONE
		INTEGER, PARAMETER :: dp=KIND(0.0D0)
		COMPLEX, PARAMETER :: j = (0, 1)
		COMPLEX(KIND=dp) :: integral
		COMPLEX(KIND=dp), DIMENSION(0:n_freq) :: f_image
		REAL(KIND=dp) :: Pi
		REAL(KIND=dp), DIMENSION(0:n_freq) :: freq
		REAL(KIND=dp), DIMENSION(0:n_time) :: correlation
		REAL(KIND=dp), DIMENSION(0:n_time) :: time
		INTEGER :: I,clock,n_freq,n_time

        !f2py intent(in) :: correlation
        !f2py intent(hide), depend(N) :: n_time = len(correlation), n_freq = len(freq)
        !f2py intent(out) :: f_image
   		Pi=4.0D0*ATAN(1.0D0)
		!omega = 2.0*Pi/n_time
		DO I=0,n_freq
          integral=(0.0_dp,0.0_dp)
          DO clock=0,n_time-1
            integral=integral+correlation(clock)*(cos(time(clock)*freq(I)*2*Pi) +j*sin(time(clock)*freq(I)*2*Pi))
          ENDDO
          f_image(I) = integral
        ENDDO
END SUBROUTINE

!def numFLap(funct, time, freq):
!    # Returns the image of $funct in $freq area.Where $funct is function of $time( in time domain).
!    # ang_freq = freq. * 2 * pi;
!    ang_freq = 2*np.pi*freq
!
!    FreqImage = np.zeros(len(ang_freq), dtype = 'complex_')
!    for omega in range(len(ang_freq)):
!        integral = 0 + 0j
!        for clock in range(len(time) - 2):
!            # preforms the Fourier - Laplace transform
!            integral += np.exp(1j * ang_freq[omega] * time[clock]) * dir_phi_n[clock]
!        FreqImage[omega] = integral
!    return FreqImage