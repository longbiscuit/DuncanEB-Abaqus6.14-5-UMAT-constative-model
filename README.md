# DuncanEB-Abaqus6.14-5-UMAT-constative-model
This is a nonlinear elastic constitutive model(Duncan-CHANG EB) UMAT widely used in the field of geotechnical engineering. It can be used in abaqus6.14-5 to calculate the dam settlement deformation.


1. 将本文件夹放在D盘，
   或修改 duncan2.py 中路径：
   os.chdir(r"D:\duncanEB")#记得将"D:\duncanEB" 引号中路径改为自己当前路径

2. 双击 main-duncan.bat 运行即可自动运行abaqus6.14-5计算，并提取CD试验仿真结果到duncanAbaquse1ev.dat和duncanAbaquse1q.dat，如果想改变加载增量步，修改loade2.dat中内容即可（这个文档在duncan.inp中 *Amplitude, name=Amp-1, INPUT=loade2.dat 使用了）

理论计算值的修正有：

    SL修正：
    if(SL>0.99)
        SL=0.99
    end
	泊松比修正：
	    v=0.5*(1-Etp/(3.0*B));
    %修正泊松比
    if(v<0.00505051)
        v=0.00505051
    end
    if(v>0.490196 )
        v=0.490196
    end
	仿真时修正除上面外，还有简单的拉裂修正，详见UMAT。


3. 实在不行，就手动操作吧。验证结果：
![duncan-e1-q](https://user-images.githubusercontent.com/21994802/192254423-59b54a90-9813-4d77-bdde-7bbe70182532.png)
![duncan-e1-ev](https://user-images.githubusercontent.com/21994802/192254464-655c0d67-ad8a-4d94-b668-9620a2c07f6c.png)


4.试验情况：
室内常规三轴排水剪切试验（CD）:
围压 3000kPa：
sc=3000;% 围压
参数：
%   K n Rf c fai dfai Kb m Kur Pa s1Ms3max s3max ssm  ，后面三个参数这个版本的umat不用了，可以不写，写了也没用。
 x=[1075,  0.33,  0.77,  0.0,  46.9,  6.5,  550,  0.25,  2150,  101.325,  0.0,  0.0,  0.0,]; %DuncanEB本构模型参数
 
 5. matlab代码计算理论值 

``` matlab
clc;
clear all;
close all;
para=[1075,  0.33,  0.77,  0.0,  46.9,  6.5,  550,  0.25,  2150,  101.325,  0.0,  0.0,  0.0];
sig30=3000;%围压
de1=ones(1,3000)*0.0001;% e1=30%，如果要进行加卸载，请给定de1合适的正负值，比如将load2中第2列各元素作差作为de1
[te1,tq,tev]=duncanModel(para,sig30,de1);
result=[te1',tq',tev'];
figure(1)
plot(result(:,1),result(:,2));
xlabel('e1');
ylabel('q');
title('EB:e1-q')
figure(2)
plot(result(:,1),result(:,3));
xlabel('e1');
ylabel('ev');
title('EB:e1-ev')


%% 子函数 进行了修正，导致后半段为直线
function  [te1,tq,tev]=duncanModel(para,sig30,de1)
%K n Rf c fai dfai Kb m Kur Pa 
 K=para(1);N=para(2);Rf=para(3);C=para(4);PHI=para(5);DPHI=para(6);Kb=para(7);m=para(8);Kur=para(9);PA=para(10);%邓肯张参数
 sc=[sig30];e0=0.75;
% e0 随便给

myType=3;%1代表不排水，2代表等p，3代表围压不变，4代表侧限,5.等向压缩

[~,number]=size(de1);%epsilon1 个数
S1=zeros(1,number);%主应力 σ1
S2=zeros(1,number);
S3=zeros(1,number);
S11=zeros(1,number);% 变换应力空间主应力 σ1
S22=zeros(1,number);
S33=zeros(1,number);
p=zeros(1,number); %变换应力空间平均应力p
q=zeros(1,number);
qz=zeros(1,number);%真实应力空间广义剪应力q
pz=zeros(1,number);
e1=zeros(1,number);%ε1
e2=zeros(1,number);
e3=zeros(1,number);
ev=zeros(1,number);%总的体变εv
evp=zeros(1,number);%塑性体积应变εv
ed=zeros(1,number);%εd 剪切应变
n=zeros(1,number);%应力比yita
nz=zeros(1,number);%真实应力空间应力比yita
e=zeros(1,number);%孔隙比
E=zeros(1,number);%弹性模量
u=zeros(1,number);%孔压
kesi=zeros(1,number);
Mf=zeros(1,number);
Mc=zeros(1,number);
Oute1=zeros(1,number);
Outq=zeros(1,number);
Outev=zeros(1,number);
PosionRatio=zeros(1,number);
SSmax=0.0;v=0.0;
b=[0];%b值，三轴b=1三轴伸长。b=0三轴压缩，b=0.25 平面应变

% sc=[400];%丰浦砂人工合成数据的围压 100 200 500 1000 
% scAlf=sprintfc('%g',sc);%把数字数组转化成字符串数组
% e0=[0.85];%丰丰浦砂人工合成数据的孔隙比 0.55 0.65 0.75 0.85
scNum=length(sc);

for im=1:scNum
    
S1(1)=sc(im)+0.000000000001;S2(1)=sc(im);S3(1)=sc(im);%三个主应力
 e(1)=e0(im);%初始孔隙比
 for in=1:number    
%     I1=S1(in)+S2(in)+S3(in);%S1真实应力主应力sig1  I1真实应力空间中第一应力不变量
%     I2=S1(in)*S2(in)+S2(in)*S3(in)+S3(in)*S1(in);
%     I3=S1(in)*S2(in)*S3(in);%没有剪应力情况下
%     F=sqrt((I1*I2-I3)/(I1*I2-9*I3));
%     qc=2*I1/(3*F-1);%  为smp准则对应的广义剪应力
%     q(in)=2*I1/(3*F-1);%变换应力空间中q，为啥和qc相同？？？？
%     p(in)=I1/3;%变换应力空间中p
%     qz(in)=1/sqrt(2)*sqrt((S1(in)-S2(in))^2+(S2(in)-S3(in))^2+(S1(in)-S3(in))^2);%真实应力空间广义剪应力
%     pz(in)=(S1(in)+S2(in)+S3(in))/3;%真实应力p，和p(in)=I1/3 一样啊！

%     S11=pz(in)+qc*(S1(in)-pz(in))/qz(in);%S11变换应力空间第一主应力，老师的经典公式
%     S22=pz(in)+qc*(S2(in)-pz(in))/qz(in);
%     S33=pz(in)+qc*(S3(in)-pz(in))/qz(in);
%     n(in)=q(in)/p(in);
%     nz(in)=qz(in)/pz(in);
    
    PHIP=PHI-DPHI*log10(S3(1)/PA);
    SL=((1-sin(PHIP/180.0*pi))*(S1(in)-S3(in)))/(2*C*cos(PHIP/180.0*pi)+2*S3(in)*sin(PHIP/180.0*pi));
    if(SL>0.99)
        SL=0.99
    end
    Et=K*PA*((S3(in)/PA)^N)*(1-Rf*SL)^2;
    Eur=Kur*PA*(S3(in)/PA)^N;
    B=Kb*PA*(S3(in)/PA)^m;
    SS=SL*((S3(in)/PA)^(0.25));

    if(SS>=SSmax)
        Etp=Et;
    elseif(SS<=0.75*SSmax)
        Etp=Eur;
    elseif(0.75*SSmax<=SS && SS<SSmax)
        Etp=Et+((SSmax-SS)/(0.25*SSmax))*(Eur-Et);
    end
    if(SSmax<SS)
    SSmax=SS;
    end

    v=0.5*(1-Etp/(3.0*B));
    %修正泊松比
    if(v<0.00505051)
        v=0.00505051
    end
    if(v>0.490196 )
        v=0.490196
    end
    PosionRatio(in)=v;
    A=1.0/Etp;
    CII=1*A;
    CIJ=-v*A;
    C11=CII;%柔度矩阵 
    C12=CIJ;
    C13=CIJ;
    C21=CIJ;
    C22=CII;
    C23=CIJ;
    C31=CIJ;
    C32=CIJ;
    C33=CII;
    %判断是哪种实验
   if myType==1
% % % % %     师姐的三轴不排水（不同的b）cu
    dS1=de1(in)/((C11+b(im)*C12)-(C11+C21+C31+b(im)*(C12+C22+C32))*((1-b(im))*C12+C13)/(C13+C23+C33+(1-b(im))*(C12+C22+C32)));
    dS3=-(C11+C21+C31+b(im)*(C12+C22+C32))*dS1/(C13+C23+C33+(1-b(im))*(C12+C22+C32));
    dS2=b(im)*dS1+(1-b(im))*dS3;
    de2=(C21+b(im)*C22)*dS1+((1-b(im))*C22+C23)*dS3;
    de3=(C31+b(im)*C32)*dS1+((1-b(im))*C32+C33)*dS3;
   elseif myType==2
        %我的等p
        dS1=(b(im)-2)*de1(in)/(-2*C11+b(im)*C11+C12-2*b(im)*C12+C13+b(im)*C13);
        dS2=(2*b(im)-1)*de1(in)/(2*C11-b(im)*C11-C12+2*b(im)*C12-C13-b(im)*C13);
        dS3=(1+b(im))*de1(in)/(-2*C11+b(im)*C11+C12-2*b(im)*C12+C13+b(im)*C13);
        de2=C21*dS1-C22*dS1-C22*dS3+C23*dS3;
        de3=C31*dS1-C32*dS1-C32*dS3+C33*dS3;
        %师姐的三轴等p
%     dS3=de1(in)/((b(im)-2)/(1+b(im))*C11+(1-2*b(im))/(1+b(im))*C12+C13);
%     dS1=(b(im)-2)/(1+b(im))*dS3;
%     dS2=(1-2*b(im))/(1+b(im))*dS3;
%     de2=C21*dS1+C22*dS2+C23*dS3;
%     de3=C31*dS1+C32*dS2+C33*dS3;
   elseif myType==3
        %我的围压不变CD
        dS1=de1(in)/(C11+b(im)*C12);
        dS2=b(im)*de1(in)/(C11+b(im)*C12);
        dS3=0;
        de2=(C21+b(im)*C22)*dS1;
        de3=(C31+b(im)*C32)*dS1;
% % %     师姐的三轴排水（不同的b）
%     dS1=de1(in)/(C11+b(im)*C12);
%     dS2=b(im)*dS1;
%     dS3=0;
%     de2=(C21+b(im)*C22)*dS1;
%     de3=(C31+b(im)*C32)*dS1;
   elseif myType==4
       %侧限压缩 
       dS1=(C22- b(im)* C22 + C23)* de1(in)/( b(im)* (C12* C21-C11* C22-C13 *C22+C12 *C23) - C13* C21 + C11* C22  + C11 *C23 -C12 *C21 );
       dS2=(-C21 + b(im)*(C21+C23)  )*de1(in)/( b(im)*(C12*C21-  C11*C22-  C13*C22+  C12*C23)  - C13*C21 + C11*C22   + C11*C23 -C12*C21 );
       dS3=(C21 + b(im)* C22)* de1(in)/( b(im)*( C12* C21+ C11*C22+ C13*C22- C12*C23) + C13*C21 - C11*C22   - C11*C23 +C12*C21);
       de2=0;
       de3=0;  
       %  师姐的 K0固结
%     dS1=de1(in)/(C11+b(im)*C12+(C13-b(im)*C12+C12)*((1-C(im))*C11+C21+C31+b(im)*((1-C(im))*C12+C22+C32))/((C(im)-1)*C13-C23-C33+(1-b(im))*((C(im)-1)*C12-C22-C32)));
%     dS3=((1-C(im))*C11+C21+C31+b(im)*((1-C(im))*C12+C22+C32))*dS1/((C(im)-1)*C13-C23-C33+(1-b(im))*((C(im)-1)*C12-C22-C32));
%     dS2=(dS1-dS3)*b(im)+dS3; 
%     de2=C21*dS1+C22*dS2+C23*dS3;
%     de3=C31*dS1+C32*dS2+C33*dS3;
      elseif myType==5
       %等向压缩 dS1==dS2==dS3
       dS1=de1(in)/(C11+C12+C13);
       dS2=dS1;
       dS3=dS1;
       de2=de1(in);
       de3=de1(in);
%        de2=(C21+C22+C23)*dS1;
%        de3=(C31+C32+C33)*dS1; 
  elseif myType==6
       %减载的三轴伸长
       dS1=de1(in)/C11;
       dS2=0;
       dS3=0;
       de2=C21*dS1;
       de3=C31*dS1;
   end 
    S1(in+1)=S1(in)+dS1;
    S2(in+1)=S2(in)+dS2;

    S3(in+1)=S3(in)+dS3;
    
    %平面应变
%     S2(in+1)=sqrt(S1(in+1)*S3(in+1));
    
    e1(in+1)=e1(in)+de1(in);
    e2(in+1)=e2(in)+de2;
    e3(in+1)=e3(in)+de3;
    ev(in+1)=e1(in+1)+e2(in+1)+e3(in+1);
    ed(in+1)=sqrt(2)/3*sqrt((e1(in+1)-e2(in+1))^2+(e1(in+1)-e3(in+1))^2+(e3(in+1)-e2(in+1))^2);
    e(in+1)=e0(im)-ev(in+1)*(1+e0(im));%孔隙比
    if(myType==1)   
%         u(in+1)=S3(1)-S3(in+1);%这个是不对的
%     qz(in)=1/sqrt(2)*sqrt((S1(in)-S2(in))^2+(S2(in)-S3(in))^2+(S1(in)-S3(in))^2);%真实应力空间广义剪应力
%     pz(in)=(S1(in)+S2(in)+S3(in))/3;%真实应力p，和p(in)=I1/3 一样啊！
     u(in+1)=(1.0/sqrt(2.0)*sqrt((S1(in+1)-S2(in+1))^2+(S2(in+1)-S3(in+1))^2+(S1(in+1)-S3(in+1))^2))/3.0+sc(im)-(S1(in+1)+S2(in+1)+S3(in+1))/3.0;
    else
        u(in+1)=0.0;
    end 
    if( in==number)
        qz(in+1) = 1.0 / sqrt(2.0) * sqrt(((S1(in+1) - S2(in+1))^2) + ((S2(in+1) - S3(in+1))^2) + ((S1(in+1) - S3(in+1))^2));%%真实应力空间广义剪应力
        pz(in+1) = (S1(in+1) + S2(in+1) + S3(in+1)) / 3.0;%%
    end 
  for un=1:number
    if(pz(un)<0)
        num=un-1;
        break;
    end
  end

 end
end
 

  te1=e1;
  tq=S1-S3;
  tev=ev;

end
```
![image](https://user-images.githubusercontent.com/21994802/192254062-0f3b7687-f649-48cf-8c87-50984d131b60.png)
![image](https://user-images.githubusercontent.com/21994802/192254123-ea305943-2c68-4eef-adbd-a5343d0595b6.png)

