function one_band_N_3_spin_symmetric_old(g0_1_1_1,g0_1_2_1,g0_1_3_1,g0_1_1_2,g0_1_2_2,g0_1_3_2,g0_1_1_3,g0_1_2_3,g0_1_3_3,u1,u2,u3,u4)
M240=g0_1_2_2^2
M120=-g0_1_2_1*g0_1_2_2
M252=M120*g0_1_1_2 + M240*g0_1_1_1
M204=-g0_1_2_1*g0_1_1_2 + g0_1_2_2*g0_1_1_1
M195=g0_1_2_2*g0_1_1_1
M60=g0_1_2_2*g0_1_1_1
M207=M204*g0_1_1_1
M51=-g0_1_2_1*g0_1_1_2 + g0_1_2_2*g0_1_1_1
M210=-g0_1_2_1*g0_1_2_2
M243=M210*g0_1_1_2 + M240*g0_1_1_1
M222=-M204*g0_1_2_1
M255=M222*g0_1_1_2 + M252*g0_1_1_1
M30=-g0_1_2_1*g0_1_1_1
M63=M30*g0_1_1_2 + M60*g0_1_1_1
M15=g0_1_1_1^2
M3264=-g0_1_3_2*g0_1_2_3 + g0_1_3_3*g0_1_2_2
M3144=g0_1_3_1*g0_1_2_3 - g0_1_3_3*g0_1_2_1
M1224=-g0_1_3_1*g0_1_2_2 + g0_1_3_2*g0_1_2_1
M3276=M1224*g0_1_1_3 + M3144*g0_1_1_2 + M3264*g0_1_1_1
M3279=M3276*g0_1_1_1
M3840=g0_1_3_3^2
M1920=-g0_1_3_2*g0_1_3_3
M4032=M1920*g0_1_2_3 + M3840*g0_1_2_2
M3552=-M3264*g0_1_3_2
M4080=M3552*g0_1_2_3 + M4032*g0_1_2_2
M3360=-g0_1_3_2*g0_1_3_3
M3888=M3360*g0_1_2_3 + M3840*g0_1_2_2
M2880=g0_1_3_3*g0_1_2_3
M2400=g0_1_3_2*g0_1_2_3
M2928=-M2400*g0_1_2_3 + M2880*g0_1_2_2
M3960=M2928*g0_1_3_1 - M3888*g0_1_2_1
M1440=-g0_1_3_2^2
M1968=-M1440*g0_1_2_3 + M1920*g0_1_2_2
M960=g0_1_3_3*g0_1_2_2
M480=-g0_1_3_2*g0_1_2_2
M1008=M480*g0_1_2_3 + M960*g0_1_2_2
M2040=-M1008*g0_1_3_1 - M1968*g0_1_2_1
M4092=M2040*g0_1_1_3 + M3960*g0_1_1_2 + M4080*g0_1_1_1
M816=-g0_1_3_2*g0_1_2_3 + g0_1_3_3*g0_1_2_2
M1848=-M816*g0_1_3_1
M3900=M1848*g0_1_1_3 + M3888*g0_1_1_1
M3912=M2880*g0_1_3_1 - M3840*g0_1_2_1
M1992=-M1920*g0_1_2_1 - M960*g0_1_3_1
M4044=M1992*g0_1_1_3 + M3912*g0_1_1_2 + M4032*g0_1_1_1
M3534=-M3276*g0_1_3_1
M4047=M3534*g0_1_1_3 + M4044*g0_1_1_1
M3120=g0_1_3_3*g0_1_2_2
M3522=-M3264*g0_1_3_1
M4035=M3522*g0_1_1_3 + M4032*g0_1_1_1
M3330=-g0_1_3_1*g0_1_3_3
M3843=M3330*g0_1_1_3 + M3840*g0_1_1_1
M771=-g0_1_3_1*g0_1_1_3 + g0_1_3_3*g0_1_1_1
M3075=g0_1_3_3*g0_1_1_1
M3084=-g0_1_3_1*g0_1_1_3 + g0_1_3_3*g0_1_1_1
M3087=M3084*g0_1_1_1
M1800=-g0_1_3_1*g0_1_3_3
M3852=M1800*g0_1_1_3 + M3840*g0_1_1_1
M3342=-M3084*g0_1_3_1
M3855=M3342*g0_1_1_3 + M3852*g0_1_1_1
M3600=g0_1_3_3*g0_1_2_3
M1560=-g0_1_3_1*g0_1_2_3
M3612=M1560*g0_1_1_3 + M3600*g0_1_1_1
M3870=M3612*g0_1_3_1 - M3852*g0_1_2_1
M1320=g0_1_3_2*g0_1_3_1
M3372=M1320*g0_1_1_3 + M3360*g0_1_1_1
M1080=-g0_1_3_1*g0_1_2_2
M3132=M1080*g0_1_1_3 + M3120*g0_1_1_1
M3390=-M3132*g0_1_3_1 - M3372*g0_1_2_1
M3903=M3390*g0_1_1_3 + M3870*g0_1_1_2 + M3900*g0_1_1_1
M3312=M3264*g0_1_2_2
M2160=g0_1_2_2*g0_1_2_3
M3192=M2160*g0_1_3_1 - M3120*g0_1_2_1
M1200=-g0_1_3_2*g0_1_2_2
M1272=-M1200*g0_1_2_1 - M240*g0_1_3_1
M3324=M1272*g0_1_1_3 + M3192*g0_1_1_2 + M3312*g0_1_1_1
M3294=-M3276*g0_1_2_1
M3327=M3294*g0_1_1_2 + M3324*g0_1_1_1
M3267=M3264*g0_1_1_1
M720=g0_1_2_2*g0_1_2_3
M978=M720*g0_1_3_1 - M960*g0_1_2_1
M498=-M240*g0_1_3_1 - M480*g0_1_2_1
M1011=M1008*g0_1_1_1 + M498*g0_1_1_3 + M978*g0_1_1_2
M828=M816*g0_1_1_1
M840=-g0_1_3_3*g0_1_2_1
M972=M840*g0_1_1_2 + M960*g0_1_1_1
M462=-M204*g0_1_3_1
M975=M462*g0_1_1_3 + M972*g0_1_1_1
M3792=M3264*g0_1_2_3
M2640=-g0_1_2_3^2
M3672=M2640*g0_1_3_1 + M3600*g0_1_2_1
M1680=g0_1_3_2*g0_1_2_3
M1752=M1680*g0_1_2_1 - M720*g0_1_3_1
M3804=M1752*g0_1_1_3 - M3672*g0_1_1_2 + M3792*g0_1_1_1
M4062=M3804*g0_1_3_1 - M4044*g0_1_2_1
M3432=M2400*g0_1_3_1 + M3360*g0_1_2_1
M1512=M1440*g0_1_2_1 - M480*g0_1_3_1
M3564=M1512*g0_1_1_3 - M3432*g0_1_1_2 + M3552*g0_1_1_1
M3582=-M3324*g0_1_3_1 - M3564*g0_1_2_1
M4095=M3582*g0_1_1_3 + M4062*g0_1_1_2 + M4092*g0_1_1_1
M786=g0_1_3_1*g0_1_2_3 - g0_1_3_3*g0_1_2_1
M306=-g0_1_3_1*g0_1_2_2 + g0_1_3_2*g0_1_2_1
M819=M306*g0_1_1_3 + M786*g0_1_1_2 + M816*g0_1_1_1
M780=g0_1_3_3*g0_1_1_1
M3090=-g0_1_3_3*g0_1_2_1
M3123=M3090*g0_1_1_2 + M3120*g0_1_1_1
M3858=M3600*g0_1_3_1 - M3840*g0_1_2_1
M3378=-M3120*g0_1_3_1 - M3360*g0_1_2_1
M3891=M3378*g0_1_1_3 + M3858*g0_1_1_2 + M3888*g0_1_1_1
M450=-g0_1_3_1*g0_1_2_2
M963=M450*g0_1_1_3 + M960*g0_1_1_1
M3282=-M3264*g0_1_2_1
M3315=M3282*g0_1_1_2 + M3312*g0_1_1_1
M888=-M816*g0_1_2_1
M1020=M1008*g0_1_1_1 + M888*g0_1_1_2
M600=g0_1_2_1*g0_1_2_3
M732=-M600*g0_1_1_2 + M720*g0_1_1_1
M990=M732*g0_1_3_1 - M972*g0_1_2_1
M360=-g0_1_3_2*g0_1_2_1
M492=-M360*g0_1_1_2 + M480*g0_1_1_1
M510=-M252*g0_1_3_1 - M492*g0_1_2_1
M1023=M1020*g0_1_1_1 + M510*g0_1_1_3 + M990*g0_1_1_2
M4050=M3792*g0_1_3_1 - M4032*g0_1_2_1
M3570=-M3312*g0_1_3_1 - M3552*g0_1_2_1
M4083=M3570*g0_1_1_3 + M4050*g0_1_1_2 + M4080*g0_1_1_1
M540=g0_1_2_3*g0_1_1_1
M798=M540*g0_1_3_1 - M780*g0_1_2_1
M300=-g0_1_3_2*g0_1_1_1
M318=-M300*g0_1_2_1 - M60*g0_1_3_1
M831=M318*g0_1_1_3 + M798*g0_1_1_2 + M828*g0_1_1_1
M270=-g0_1_3_1*g0_1_1_1
M783=M270*g0_1_1_3 + M780*g0_1_1_1
M3102=-M3084*g0_1_2_1
M3135=M3102*g0_1_1_2 + M3132*g0_1_1_1
M237=M204*g0_1_1_2
M225=g0_1_2_2*g0_1_1_2
M45=g0_1_1_1*g0_1_1_2
M561=g0_1_2_2*g0_1_1_3 - g0_1_2_3*g0_1_1_2
M765=M252*g0_1_1_3 - M732*g0_1_1_2
M717=M204*g0_1_1_3
M573=-M540*g0_1_1_2 + M60*g0_1_1_3
M753=M240*g0_1_1_3 - M720*g0_1_1_2
M525=g0_1_1_1*g0_1_1_3
M705=g0_1_2_2*g0_1_1_3
M723=M210*g0_1_1_3 + M720*g0_1_1_1
M531=-g0_1_2_1*g0_1_1_3 + g0_1_2_3*g0_1_1_1
M735=M222*g0_1_1_3 + M732*g0_1_1_1
M543=M30*g0_1_1_3 + M540*g0_1_1_1
M291=g0_1_3_1*g0_1_1_2 - g0_1_3_2*g0_1_1_1
M483=-M450*g0_1_1_2 + M480*g0_1_1_1
M495=-M462*g0_1_1_2 + M492*g0_1_1_1
M303=-M270*g0_1_1_2 + M300*g0_1_1_1
nn2=-u1^2*M3840 + u1^2*g0_1_3_3 - (-u1 + u2)^2*M3891 + (-u1 + u2)^2*M819 - (-u1 + u3)^2*M4044 + (-u1 + u3)^2*M972 + (u1 - u2 - u3 + u4)^2*M1023 - (u1 - u2 - u3 + u4)^2*M4095 - (-u1 + u2)*u1*M3843 - (-u1 + u2)*u1*M3888 + (-u1 + u2)*u1*M771 + (-u1 + u2)*u1*M816 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M1011 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M3903 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M4083 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M831 - (-u1 + u3)*u1*M3852 - (-u1 + u3)*u1*M4032 + (-u1 + u3)*u1*M780 + (-u1 + u3)*u1*M960 - (-u1 + u3)*(-u1 + u2)*M3900 - (-u1 + u3)*(-u1 + u2)*M4035 + (-u1 + u3)*(-u1 + u2)*M828 + (-u1 + u3)*(-u1 + u2)*M963 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M1020 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M4047 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M4092 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M975 + (u1 - u2 - u3 + u4)*u1*M1008 - (u1 - u2 - u3 + u4)*u1*M3855 - (u1 - u2 - u3 + u4)*u1*M4080 + (u1 - u2 - u3 + u4)*u1*M783
g_1_3=M561*((-u1 + u2)*u1 + (-u1 + u2)^2) + M705*((-u1 + u3)*u1 + (-u1 + u3)*(-u1 + u2)) + M717*((-u1 + u3)*(u1 - u2 - u3 + u4) + (-u1 + u3)^2) + M753*((-u1 + u2)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)*u1) + M765*((-u1 + u3)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)^2) + ((-u1 + u2)*(u1 - u2 - u3 + u4) + (-u1 + u3)*(-u1 + u2))*M573 + ((-u1 + u3)*u1 + (u1 - u2 - u3 + u4)*u1)*M525 + ((-u1 + u2)*u1 + u1^2)*g0_1_1_3
g_2_1=M222*(-(-u1 + u3)*(u1 - u2 - u3 + u4) - (-u1 + u3)^2) + M30*(-(-u1 + u3)*u1 - (-u1 + u3)*(-u1 + u2)) + (-(-u1 + u3)*u1 - (u1 - u2 - u3 + u4)*u1)*M210 - (-(-u1 + u2)*u1 - u1^2)*g0_1_2_1
nn1=u1^2*M3840 - 2*u1^2*g0_1_3_3 - (-u1 + u2)^2*M3123 + (-u1 + u2)^2*M3891 + (-u1 + u2)^2*M51 - (-u1 + u2)^2*M819 + (-u1 + u3)^2*M204 - (-u1 + u3)^2*M3276 + (-u1 + u3)^2*M4044 - (-u1 + u3)^2*M972 - (u1 - u2 - u3 + u4)^2*M1023 + (u1 - u2 - u3 + u4)^2*M255 - (u1 - u2 - u3 + u4)^2*M3327 + (u1 - u2 - u3 + u4)^2*M4095 - (-u1 + u2)*u1*M3075 - (-u1 + u2)*u1*M3120 + (-u1 + u2)*u1*M3843 + (-u1 + u2)*u1*M3888 - (-u1 + u2)*u1*M771 - (-u1 + u2)*u1*M816 + (-u1 + u2)*u1*g0_1_1_1 + (-u1 + u2)*u1*g0_1_2_2 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M1011 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M243 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M3135 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M3315 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M3903 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M4083 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M63 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M831 - (-u1 + u3)*u1*M3084 - (-u1 + u3)*u1*M3264 + (-u1 + u3)*u1*M3852 + (-u1 + u3)*u1*M4032 - (-u1 + u3)*u1*M780 - (-u1 + u3)*u1*M960 + (-u1 + u3)*u1*g0_1_1_1 + (-u1 + u3)*u1*g0_1_2_2 + (-u1 + u3)*(-u1 + u2)*M195 - (-u1 + u3)*(-u1 + u2)*M3132 - (-u1 + u3)*(-u1 + u2)*M3267 + (-u1 + u3)*(-u1 + u2)*M3900 + (-u1 + u3)*(-u1 + u2)*M4035 + (-u1 + u3)*(-u1 + u2)*M60 - (-u1 + u3)*(-u1 + u2)*M828 - (-u1 + u3)*(-u1 + u2)*M963 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M1020 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M207 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M252 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M3279 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M3324 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M4047 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M4092 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M975 - (u1 - u2 - u3 + u4)*u1*M1008 + (u1 - u2 - u3 + u4)*u1*M15 + (u1 - u2 - u3 + u4)*u1*M240 - (u1 - u2 - u3 + u4)*u1*M3087 - (u1 - u2 - u3 + u4)*u1*M3312 + (u1 - u2 - u3 + u4)*u1*M3855 + (u1 - u2 - u3 + u4)*u1*M4080 - (u1 - u2 - u3 + u4)*u1*M783 + u1^2
nn3=-u1^2*M3840 + u1^2*g0_1_3_3 + (-u1 + u2)^2*M3123 - (-u1 + u2)^2*M3891 + (-u1 + u3)^2*M3276 - (-u1 + u3)^2*M4044 + (u1 - u2 - u3 + u4)^2*M3327 - (u1 - u2 - u3 + u4)^2*M4095 + (-u1 + u2)*u1*M3075 + (-u1 + u2)*u1*M3120 - (-u1 + u2)*u1*M3843 - (-u1 + u2)*u1*M3888 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M3135 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M3315 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M3903 - (-u1 + u2)*(u1 - u2 - u3 + u4)*M4083 + (-u1 + u3)*u1*M3084 + (-u1 + u3)*u1*M3264 - (-u1 + u3)*u1*M3852 - (-u1 + u3)*u1*M4032 + (-u1 + u3)*(-u1 + u2)*M3132 + (-u1 + u3)*(-u1 + u2)*M3267 - (-u1 + u3)*(-u1 + u2)*M3900 - (-u1 + u3)*(-u1 + u2)*M4035 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M3279 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M3324 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M4047 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M4092 + (u1 - u2 - u3 + u4)*u1*M3087 + (u1 - u2 - u3 + u4)*u1*M3312 - (u1 - u2 - u3 + u4)*u1*M3855 - (u1 - u2 - u3 + u4)*u1*M4080
nn4=u1^2*M3840 + (-u1 + u2)^2*M3891 + (-u1 + u3)^2*M4044 + (u1 - u2 - u3 + u4)^2*M4095 + (-u1 + u2)*u1*M3843 + (-u1 + u2)*u1*M3888 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M3903 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M4083 + (-u1 + u3)*u1*M3852 + (-u1 + u3)*u1*M4032 + (-u1 + u3)*(-u1 + u2)*M3900 + (-u1 + u3)*(-u1 + u2)*M4035 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M4047 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M4092 + (u1 - u2 - u3 + u4)*u1*M3855 + (u1 - u2 - u3 + u4)*u1*M4080
g_2_2=M252*((-u1 + u3)*(u1 - u2 - u3 + u4) + (-u1 + u3)^2) + M255*((-u1 + u3)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)^2) + M51*((-u1 + u2)*u1 + (-u1 + u2)^2) + M60*((-u1 + u3)*u1 + (-u1 + u3)*(-u1 + u2)) + M63*((-u1 + u2)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)*u1) + ((-u1 + u2)*(u1 - u2 - u3 + u4) + (-u1 + u3)*(-u1 + u2))*M243 + ((-u1 + u3)*u1 + (u1 - u2 - u3 + u4)*u1)*M240 + ((-u1 + u2)*u1 + u1^2)*g0_1_2_2
g_3_2=u1^2*g0_1_3_2 - (-u1 + u3)^2*M492 - (-u1 + u2)*u1*M291 - (-u1 + u3)*u1*M300 - (-u1 + u3)*u1*M480 - (-u1 + u3)*(-u1 + u2)*M483 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M495 - (u1 - u2 - u3 + u4)*u1*M303
g_3_1=u1^2*g0_1_3_1 - (-u1 + u3)^2*M462 - (-u1 + u2)*u1*M306 - (-u1 + u3)*u1*M270 - (-u1 + u3)*u1*M450 - (-u1 + u3)*(-u1 + u2)*M318 - (-u1 + u3)*(u1 - u2 - u3 + u4)*M510 - (u1 - u2 - u3 + u4)*u1*M498
g_2_3=M531*((-u1 + u2)*u1 + (-u1 + u2)^2) + M540*((-u1 + u3)*u1 + (-u1 + u3)*(-u1 + u2)) + M543*((-u1 + u2)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)*u1) + M732*((-u1 + u3)*(u1 - u2 - u3 + u4) + (-u1 + u3)^2) + M735*((-u1 + u3)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)^2) + ((-u1 + u2)*(u1 - u2 - u3 + u4) + (-u1 + u3)*(-u1 + u2))*M723 + ((-u1 + u3)*u1 + (u1 - u2 - u3 + u4)*u1)*M720 + ((-u1 + u2)*u1 + u1^2)*g0_1_2_3
g_1_1=M195*((-u1 + u3)*u1 + (-u1 + u3)*(-u1 + u2)) + M207*((-u1 + u3)*(u1 - u2 - u3 + u4) + (-u1 + u3)^2) + M243*((-u1 + u2)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)*u1) + M255*((-u1 + u3)*(u1 - u2 - u3 + u4) + (u1 - u2 - u3 + u4)^2) + M51*((-u1 + u2)*u1 + (-u1 + u2)^2) + ((-u1 + u2)*(u1 - u2 - u3 + u4) + (-u1 + u3)*(-u1 + u2))*M63 + ((-u1 + u3)*u1 + (u1 - u2 - u3 + u4)*u1)*M15 + ((-u1 + u2)*u1 + u1^2)*g0_1_1_1
g_1_2=M225*((-u1 + u3)*u1 + (-u1 + u3)*(-u1 + u2)) + M237*((-u1 + u3)*(u1 - u2 - u3 + u4) + (-u1 + u3)^2) + ((-u1 + u3)*u1 + (u1 - u2 - u3 + u4)*u1)*M45 + ((-u1 + u2)*u1 + u1^2)*g0_1_1_2
g_3_3=u1^2*g0_1_3_3 + (-u1 + u2)^2*M819 + (-u1 + u3)^2*M972 + (u1 - u2 - u3 + u4)^2*M1023 + (-u1 + u2)*u1*M771 + (-u1 + u2)*u1*M816 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M1011 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M831 + (-u1 + u3)*u1*M780 + (-u1 + u3)*u1*M960 + (-u1 + u3)*(-u1 + u2)*M828 + (-u1 + u3)*(-u1 + u2)*M963 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M1020 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M975 + (u1 - u2 - u3 + u4)*u1*M1008 + (u1 - u2 - u3 + u4)*u1*M783
p=(-u1 + u2)^2*M51 + (-u1 + u3)^2*M204 + (u1 - u2 - u3 + u4)^2*M255 + (-u1 + u2)*u1*g0_1_1_1 + (-u1 + u2)*u1*g0_1_2_2 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M243 + (-u1 + u2)*(u1 - u2 - u3 + u4)*M63 + (-u1 + u3)*u1*g0_1_1_1 + (-u1 + u3)*u1*g0_1_2_2 + (-u1 + u3)*(-u1 + u2)*M195 + (-u1 + u3)*(-u1 + u2)*M60 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M207 + (-u1 + u3)*(u1 - u2 - u3 + u4)*M252 + (u1 - u2 - u3 + u4)*u1*M15 + (u1 - u2 - u3 + u4)*u1*M240 + u1^2
[p,nn1,nn2,nn3,nn4,g_1_1,g_2_1,g_3_1,g_1_2,g_2_2,g_3_2,g_1_3,g_2_3,g_3_3]
end
