      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
CCC*******************************************************************************************
C    ————————————————————————————————————————————————————————————————————————————————————————————————
C     Author:Zhu Binglong  email:blzhu@buaa.edu.cn  
C     Solution Technique:QUASI-NEWTON 
C     2021年3月16日
C    ————————————————————————————————————————————————————————————————————————————————————————————————
CCC  自定义数组的定义       
      DIMENSION STRE(6),STRA(6),DSTRA(6),PRES(6),AN(3,3)
      DIMENSION PS(3),DSTRESS(NTENS)
      !下面的带码debug时作暂停功能使用 https://blog.csdn.net/sinat_33528967/article/details/52002321    
        !    logical :: firstrun = .true.
        !integer::tempvar,temKSTEP,temKINC,temNPT
        !if(firstrun) then
        !write(*,*)"please input an integer:"
        !read(*,*)tempvar
        !firstrun = .false.
        !endif
        !temKSTEP=KSTEP
        !temKINC =KINC
        !temNPT=NPT !设置断点 积分点
        !IF(KINC.GT.270)THEN
        ! tempvar=284
        !END IF  
CCC*******************************************************************************************
C     PNEWDT 可用来控制时间步的变化。如果设置为小于 1 的数，则程序放弃当前计算，并用新的时间增量 DTIME X PNEWDT
C     作为新的时间增量计算；这对时间相关的材料如聚合物等有用；
C      如果设为大余 1 的数，则下一个增量步加大 DTIME 为 DTIME X PNEWDT。可以更新。
CCC   读取自定义邓肯张模型本构参数，参数共有10个，在自定义材料表格中填10个参数即可  
      EK=PROPS(1)!K
      EN=PROPS(2)!n
      RF=PROPS(3)!Rf
      C=PROPS(4)+1.0D0!粘聚力c/kPa ; +5.0kPa会加快计算速度，不加的话计算超级慢。这5kPa相对于大坝中部几兆帕的围压影响很小
      FAIINI=PROPS(5)/180.0*3.141592653!内摩擦角，填入角度。会在这儿转化为弧度
      DFAI=PROPS(6)/180.0*3.141592653!剪涨角
      EKB=PROPS(7)!Kb
      EM=PROPS(8)!m
      EKUR=PROPS(9)!Kur
      PA=PROPS(10)!标准大气压 101.325kPa
      
C     单元初始应力状态，若单元第一次参加计算，则sigma3=0，此时需要强制给定一个较小的围压，可取本层土中部高度处单元竖向应力，一般可令sigma3=10kPa
      Sigma3Ini=7.25D0
      
CCC   下面是状态变量，用来记录此积分点历史上某值,至少给2个，正常给12个
      
C     1.状态变量 SSmax:历史上最大的应力状态函数值 stress state max value
      IF(STATEV(1).NE.0.0)THEN
      SSmax=STATEV(1)
      ELSE
        SSmax=0.0D0
      END IF
C     2.积分点历史上最大的围压sigma3
      IF(STATEV(2).NE.0.0)THEN 
       S3O=STATEV(2)
       ELSE
       S3O=Sigma3Ini
       END IF
       

C    另外还有密度，比如为2.3t/m3；重力加速度为-9.8m/s2；应力kPa；长度为m；单位一定要统一。一般不用缩减积分单元
CCC*******************************************************************************************
CCC  将Abaqus中弹性力学习惯的应力转化为土力学习惯的应力（土力学中压应力居多，所以以压为正，以拉为负，而弹性力学中符号正好相反）    
      DO 10 K1=1,6
	STRE(K1)=-STRESS(K1)
	STRA(K1)=-STRAN(K1)
	DSTRA(K1)=-DSTRAN(K1)         
10    CONTINUE
      
      
CCC   求当前应力状态的刚度矩阵，并根据应变增量更新应力   
      
C     求解主应力PS(1)、PS(2)、PS(3)分别代表最大主应力sigma1，中主应力sigma2和最小主应力sigma3
      CALL PRINC_STRE(STRE,PS)
C      修正主应力 PS        
      CALL MODIFYPS(KSTEP,KINC,Sigma3Ini,PS)

C    根据参数和当前sigma3求模量B和EMOD 
      CALL GETMOD(FAIINI,DFAI,PS,PA,C,SSMAX,EKB,EM,EK,EN,RF,EKUR,
     1SS,SL,S3O,B,EMOD)

C     根据B和EMOD 组装当前应力状态的弹性刚度矩阵
      CALL STIFFNESS(B,EMOD,DDSDDE)
      
C    更新应力状态
      DO 30 K1=1,6
      DO 20 K2=1,6
      STRE(K2)=STRE(K2)+DDSDDE(K2,K1)*DSTRA(K1)
20    CONTINUE
30    CONTINUE       
      
      
CCC   求更新应力状态的弹性刚度矩阵 
      
C     求更新应力的主应力      
       CALL PRINC_STRE(STRE,PS)
C      修正主应力 PS       
       CALL MODIFYPS(KSTEP,KINC,Sigma3Ini,PS)
      
C    根据参数和当前sigma3求模量B和EMOD    
      CALL GETMOD(FAIINI,DFAI,PS,PA,C,SSMAX,EKB,EM,EK,EN,RF,EKUR,
     1SS,SL,S3O,B,EMOD)
      
C     求更新应力的弹性刚度矩阵
      CALL STIFFNESS(B,EMOD,DDSDDE)
      
      
CCC   更新状态变量      
      IF(SS.GT.SSmax)SSmax=SS
      IF(PS(3).GT.S3O)S3O=PS(3)
      
      
CCC   将符号改回来，使之符合弹性力学标准（ABAQUS默认弹性力学标准）
      DO 40 K1=1,6
      STRESS(K1)=-STRE(K1)
      STRAN(K1)=-STRA(K1)
      DSTRAN(K1)=-DSTRA(K1)
40    CONTINUE
      
      
CCC   传递状态变量到本积分点下一迭代步
C     下面的状态变量必不可少：       
c     1.历史上最大的应力状态     
      STATEV(1)=SSmax
c     2.历史上最大的围压sigma3 
      STATEV(2)=S3O 
      
C     下面的状态变量可有可无，主要用来记录积分点的变量值：
c     3.应力水平SL    
      STATEV(3)=SL
c     4.EMOD
      STATEV(4)=EMOD
c     5.B         
      STATEV(5)=B
c     6.泊松比 ENU        
      STATEV(6)=0.5D0*(1.0-EMOD/3.0D0/B)
c     7. epsilonx epsilony epsilonz epsilonxy epsilonxz epsilonyz
      DO 50 K1=1,6
      STATEV(K1+6)=STRA(K1)
50    CONTINUE      

      END
 
CCC********************************************************************************************C
C======================================子例行函数===============================================C
CCC********************************************************************************************C
 
CCC*********************************************************************************************
CCC 1、计算土力学中的主应力和主应力方向
      SUBROUTINE PRINC_STRE(STRE,PS)
CCC---------------------------------------------------------------------------------------------
CCC   IN    STRE(6)——应力张量
CCC   OUT   PS(3)——主应力，PS(1)=sigma1,PS(3)=sigma3
CCC   OUT   AN(3,3)——主应力方向
CCC---------------------------------------------------------------------------------------------       
      INCLUDE 'ABA_PARAM.INC'
!      
      DIMENSION STRE(6),PS(3),AN(3,3)
!      
      CALL SPRIND(STRE,PS,AN,1,3,3)
      IF(PS(1).LT.PS(2))THEN
      TR=PS(1)
      PS(1)=PS(2)
      PS(2)=TR     
      ENDIF
      IF(PS(1).LT.PS(3))THEN
      TR=PS(1)
      PS(1)=PS(3)
      PS(3)=TR      
      ENDIF
      IF(PS(2).LT.PS(3))THEN
      TR=PS(2)
      PS(2)=PS(3)
      PS(3)=TR     
      ENDIF
      END         
CCC*********************************************************************************************
CCC 2、计算模量
      SUBROUTINE GETMOD(FAIINI,DFAI,PS,PA,C,SSMAX,EKB,EM,EK,EN,RF,EKUR,
     1SS,SL,S3O,B,EMOD)
CCC---------------------------------------------------------------------------------------------
CCC   IN    FAIINI——初始内摩擦角fai0
CCC   IN    DFAI——delta fai     
CCC   IN    PS(3)——主应力数组  
CCC   IN    PA——标准大气压力，一般为101.325kPa，消除单位用
CCC   IN    C——粘聚力
CCC   IN    SSMAX——历史上最大的应力状态函数  
CCC   IN    EKB—参数Kb
CCC   IN    EM—参数m
CCC   IN    EK—参数K
CCC   IN    EN—参数n
CCC   IN    RF—参数Rf
CCC   IN    EKUR—参数Kur
CCC   IN    SS—应力状态函数 stress state function
CCC   IN    SL—应力水平 stress level
CCC   IN    S30—历史上最大sigma3（最小主应力）
CCC   OUT   B——体积模量B  
CCC   OUT   EMOD——弹性模量Et
CCC---------------------------------------------------------------------------------------------  
      INCLUDE 'ABA_PARAM.INC'
C     声明数组      
      DIMENSION PS(3)
C     根据当前围压计算非线性摩擦角  
      FAI=FAIINI-DFAI*LOG10(PS(3)/PA)
C     判断加卸载
      
c     应力水平   
      SL=((1-SIN(FAI))*(PS(1)-PS(3)))/
     1(2.0D0*C*COS(FAI)+2.0D0*PS(3)*SIN(FAI))
      
c     SL修正，SL一般不大于0.95      
      IF(SL.GT.0.99D0)THEN
          SL=0.99D0
c          WRITE(7,*)'SL.GT.0.95'!MESSAGE
      END IF
      
c     应力状态函数SS
      SS=SL*(S3O/PA)**0.25D0
      
c     求体积模量B
      B=EKB*PA*((PS(3)/PA)**EM)
c     求加卸载模量Et和Etur
      Et=EK*PA*((PS(3)/PA)**EN)*((1.0D0-RF*SL)**2.0D0)
      Etur=EKUR*PA*((PS(3)/PA)**EN)


C     根据加卸载判断决定所使用的的弹性模量  
      IF(SS.GE.SSmax )THEN
       EMOD=Et
      ELSEIF(SS.GT.(0.75*SSmax).AND. SS.LT.SSmax)THEN
       EMOD=Et+(SSmax-SS)*(Etur-Et)/(0.25*SSmax)
      ELSEIF(SS.LE.(0.75*SSmax))THEN 
       EMOD=Etur
      END IF
     
      
C     修正模量      
      IF ( isnan(EMOD) ) THEN
        EMOD= 0.25*EK*PA*((0.02)**EN)
c       WRITE(7,*)'EMOD is NAN'! 用不着
      END IF
      
c     仅用作观察泊松比v      
c     ENU=0.5D0*(1.0-EMOD/3.0D0/B) 
c     IF(ENU.GT.0.490196)then
c      WRITE(7,*)'2-ENU.GT.0.490196'!MESSAGE
c     END IF
c     IF(ENU.LT.0.00505051)THEN
c     WRITE(7,*)'2-ENU.LT.0.00505051'!MESSAGE
c     END IF
      
c     修正B
c     下限修正：epsilon1较小时用几次，epsilon1稍微变大时就用不着了   
      IF(B .LT. (0.33*EMOD))THEN
       B=0.33D0*EMOD!保证v>0.00505051
c       WRITE(7,*)'B .LT. (0.33*EMOD)'!MESSAGE
      END IF
c     上限修正能用着，一般取0.49      
      IF(B .GT. (17.0*EMOD))THEN
       B=17.0D0*EMOD!保证v<0.490196
c       WRITE(7,*)'B .GT. (17.0*EMOD)'!MESSAGE
      END IF
      
      END

CCC*********************************************************************************************
CCC 3、计算弹性刚度矩阵
      SUBROUTINE STIFFNESS(B,Etq,DDSDDE)
CCC---------------------------------------------------------------------------------------------
CCC   IN    B——体积模量
CCC   IN    Etq——杨氏模量E      
CCC   OUT   DDSDDE(6,6)——弹性刚度矩阵[De]
CCC---------------------------------------------------------------------------------------------        
      INCLUDE 'ABA_PARAM.INC'
!      
      DIMENSION DDSDDE(6,6)

      A_TBaddEtq=(3*B/(9*B-Etq))*(3.0D0*B+Etq)
      B_TBminusEtq=(3*B/(9*B-Etq))*(3.0D0*B-Etq)
      C_Etq=(3*B/(9*B-Etq))*Etq
      DO 20 K1=1,6
      DO 10 K2=1,6
      DDSDDE(K2,K1)=0.0
10    CONTINUE
20    CONTINUE
      DO 40 K1=1,3
      DO 30 K2=1,3
      DDSDDE(K2,K1)=B_TBminusEtq
30    CONTINUE
      DDSDDE(K1,K1)=A_TBaddEtq
40    CONTINUE
      DO 50 K1=4,6
      DDSDDE(K1,K1)=C_Etq
50    CONTINUE
      END   
CCC*********************************************************************************************
CCC 4、计算弹性刚度矩阵
CCC   IN    KSTEP——当前计算步
CCC   IN    KINC——当前增量步  
CCC   IN    Sigma3Ini——新添加计算单元初始应力中的sigma3
CCC   IN/OUT   PS(3)——修正前、后的主应力  
CCC---------------------------------------------------------------------------------------------  
      SUBROUTINE MODIFYPS(KSTEP,KINC,Sigma3Ini,PS)
      INCLUDE 'ABA_PARAM.INC'
C     声明数组      
      DIMENSION PS(3)
C        首先判断是否受拉 
          IF(KINC>1 .AND. PS(3).LT.1.0D0)THEN
              !WRITE(7,*)KSTEP,KINC,'1- KINC>1, PS(3).LT.1.0E-01'!DAT
              PS(3)=1.0D0
              IF(PS(1).LT.1.0D0)THEN
              PS(1)=2.0D0*PS(3)
              PS(2)=PS(3)
              END IF
          ELSE IF(KINC.EQ.1 .AND. PS(3).LE.0.0)THEN
              !WRITE(7,*)KSTEP,KINC,'2- KINC<=1, PS(3).LT.0.0'!DAT
              PS(3)=Sigma3Ini
              IF(PS(1).LE.PS(3))THEN
              PS(1)=2.0D0*Sigma3Ini
              PS(2)=Sigma3Ini
              END IF
          END IF
      
      END
