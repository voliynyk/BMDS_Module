      DOUBLE PRECISION FUNCTION SLOGF(X)
      DOUBLE PRECISION X, COEFS(4), V, DLOG
      INTEGER I
      DATA COEFS/6.7165863851209542D+50,
     * -2.0154759155362862e+35,
     *  2.0169759155362859e+19,
     * -710/
      IF (X .GE. 1D-16) THEN
         SLOGF = DLOG(X)
      ELSE
         V = 0.0
         DO I = 1, 4
            V = X * V + COEFS(I)
         ENDDO
         SLOGF = V
      ENDIF
      RETURN
      END
