C 1 2 3 4 5 6
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C COMPUTATIONAL IMPLEMENTATION OF GAUSS ELIMINATION METHOD EXTRACTED 
C FORM RAO'S 2ND EDITION OF "THE FINITE ELEMENT METHOD IN ENGINEERING"
C
C=======================================================================
C
C     MAIN PROGRAM WHICH CALLS THE SUBROUTINE GAUSS
C
C=======================================================================
      PROGRAM MAIN
      
        DIMENSION A(3,3),B(3,1),LP(3),LQ(3,2),R(3)
     
        DATA ((A(I,J),J=1,3), I=1,3)/0.0,2.0,1.0,1.0,-2.0,-3.0,-1.0,1,2/
        DATA (B(I,1),I=1,3)/-8.0,0,3/
        DATA N,M,IFLAG/3,1,2/
     
        PRINT 10, ((A(I,J),J=1,3),I=1,3)
        PRINT 20, (B(I,1),I=1,3)
     
        CALL GAUSS(A,B,N,M,IFLAG,LP,LQ,R)
     
        PRINT 30, ((A(I,J),J=1,3),I=1,3)
        PRINT 40, (B(I,1),I=1,3)
     
 10     FORMAT (2X,'ORIGINAL COEFFICIENT MATRIX',//,3(E12.4, 1X))
 20     FORMAT (//,2X,'RIGHT HAND SIDE VECTOR',//,3(E12.4, 1X))
 30     FORMAT (//,2X,'INVERSE OF COEFFICIENT MATRIX',//,3(E12.4, 1X))
 40     FORMAT (//,2X,'SOLUTION VECTOR',//,3(E12.4, 1X))  
   
        STOP
      END PROGRAM MAIN
      
      SUBROUTINE GAUSS(A,B,N,M,IFLAG,LP,LQ,R)
        DIMENSION A(N,N),B(N,M),LP(N),LQ(N,2),R(N)
        DO 10 I=1,N
 10       LP(I)=0
        DO 150 K=1,N
          CON = 0.0
          DO 50 I=1,N
            IF (LP(I).EQ.1) GOTO 50
            DO 40 J=1,N
              IF (LP(J) - 1) 30,40,200
 30           IF (ABS(CON).GE.ABS(A(I,J))) GOTO 40
              IR=I
              IC=J
              CON=A(I,J)   
 40        CONTINUE
 50       CONTINUE
          LP(IC)=LP(IC)+1
          IF (IR.EQ.IC) GOTO 90
          DO 60 I=1,N
            CON=A(IR,I)
            A(IR,I)=A(IC,I)
 60         A(IC,I)=CON
          IF (IFLAG.EQ.0) GOTO 90
          DO 70 I=1,M
            CON=B(IR,1)
            B(IR,I)=B(IC,I)
 70         B(IC,I)=CON
 90       LQ(K,1)=IR
          LQ(K,2)=IC
          R(K)=A(IC,IC)
          A(IC,IC)=1.0
          DO 100 I=1,N
 100        A(IC,I)=A(IC,I)/R(K)
          IF (IFLAG.EQ.0) GOTO 120
          DO 110 I=1,M
 110        B(IC,1)=B(IC,1)/R(K)
 120      DO 150 I=1,N
            IF (I.EQ.IC) GOTO 150
            CON=A(I,IC)
            A(I,IC)=0.0
            DO 130 J=1,N
 130          A(I,J)=A(I,J)-A(IC,J)*CON
          IF (IFLAG.EQ.0) GOTO 150
            DO 140 II=1,M
 140          B(I,II)=B(I,II)-B(IC,II)*CON
 150      CONTINUE
       DO 170 I=1,N
         J=N-I+1
         IF (LQ(J,1).EQ.LQ(J,2)) GOTO 170
         IR=LQ(J,1)
         IC=LQ(J,2)
         DO 160 K=1,N
           CON=A(K,IC)
           A(K,IR)=A(K,IC)
           A(K,IC)=CON
 160     CONTINUE
 170    CONTINUE 
 200    RETURN
      END SUBROUTINE GAUSS
