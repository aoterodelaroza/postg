! Copyright (c) 2013 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider <hs7@post.queensu.ca>,
! and Axel D. Becke <axel.becke@dal.ca>
!
! postg is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
module param

  ! units
  integer, parameter :: iin = 5
  integer, parameter :: iout = 6
  integer, parameter :: ihrsh = 9
  integer, parameter :: luwfn = 10
  integer, parameter :: imosa = 11

  ! constants
  real*8, parameter :: pi = 3.141592653589793d0
  real*8, parameter :: fourpi = 4d0*pi
  real*8, parameter :: third = 1d0/3d0
  real*8, parameter :: third2 = 2d0/3d0
  real*8, parameter :: small = 1d-10
  integer, parameter :: mline = 1024
  
  ! chars
  character*(1), parameter :: blank = " " !< blank
  character*(1), parameter :: null = char(0) !< null character (ascii 0)
  character*(1), parameter :: newline = char(10) !< newline char (ascii 10 in unix)

  ! atomic symbols
  character*2, parameter :: ptable(0:103) = (/&
     'Bq',&
     'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',&
     'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
     'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
     'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
     'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',&
     'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',&
     'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
     'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',&
     'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',&
     'Md','No','Lr'/)
  
  ! Free Atomic Polarizabilities from
  ! CRC Handbook of Chemistry and Physics, 88th Ed.
  ! untis are AA^3, then transformed to au (0.529...)^3.
  real*8, parameter ::  frepol(0:102) = (/&
      0.d0,&
      0.6668D0,   0.2051D0,  24.3300D0,   5.6000D0,   3.0300D0,   1.7600D0,&
      1.1000D0,   0.8020D0,   0.5570D0,   0.3956D0,  24.1100D0,  10.6000D0,&
      6.8000D0,   5.3800D0,   3.6300D0,   2.9000D0,   2.1800D0,   1.6411D0,&
     43.4000D0,  22.8000D0,  17.8000D0,  14.6000D0,  12.4000D0,  11.6000D0,&
      9.4000D0,   8.4000D0,   7.5000D0,   6.8000D0,   6.2000D0,   5.7500D0,&
      8.1200D0,   6.0700D0,   4.3100D0,   3.7700D0,   3.0500D0,   2.4844D0,&
     47.3000D0,  27.6000D0,  22.7000D0,  17.9000D0,  15.7000D0,  12.8000D0,&
     11.4000D0,   9.6000D0,   8.6000D0,   4.8000D0,   7.2000D0,   7.3600D0,&
     10.2000D0,   7.7000D0,   6.6000D0,   5.5000D0,   5.3500D0,   4.0440D0,&
     59.4200D0,  39.7000D0,  31.1000D0,  29.6000D0,  28.2000D0,  31.4000D0,&
     30.1000D0,  28.8000D0,  27.7000D0,  23.5000D0,  25.5000D0,  24.5000D0,&
     23.6000D0,  22.7000D0,  21.8000D0,  21.0000D0,  21.9000D0,  16.2000D0,&
     13.1000D0,  11.1000D0,   9.7000D0,   8.5000D0,   7.6000D0,   6.5000D0,&
      5.8000D0,   5.0200D0,   7.6000D0,   6.8000D0,   7.4000D0,   6.8000D0,&
      6.0000D0,   5.3000D0,  48.6000D0,  38.3000D0,  32.1000D0,  32.1000D0,&
     25.4000D0,  24.9000D0,  24.8000D0,  24.5000D0,  23.3000D0,  23.0000D0,&
     22.7000D0,  20.5000D0,  19.7000D0,  23.8000D0,  18.2000D0,  17.5000D0/)&
     / (0.52917720859d0)**3

  !! data types
  ! molecular mesh
  type tmesh
     integer :: n
     real*8, allocatable :: w(:), x(:,:)
     real*8, allocatable :: rho(:,:), b(:,:)
  end type tmesh
  ! molecule
  type molecule
     character(mline) :: name
     integer :: n
     real*8, allocatable :: x(:,:)
     integer, allocatable :: z(:)
     integer :: nmo, npri
     integer, allocatable :: icenter(:)
     integer, allocatable :: itype(:)
     real*8, allocatable :: e(:)
     real*8, allocatable :: occ(:)
     real*8, allocatable :: c(:,:), mm(:,:), v(:), q(:)
     integer :: mult
     integer :: wfntyp  ! 0 - closed, 1 - open, 2 - restricted open, 3 - fractional
     real*8 :: nelec, charge, egauss
     logical :: useecp
  end type molecule

  ! % of HF
  real*8 :: chf 
  integer, parameter :: chf_blyp = -1
  integer, parameter :: chf_b3lyp = -2
  integer, parameter :: chf_bhahlyp = -3
  integer, parameter :: chf_camb3lyp = -4
  integer, parameter :: chf_pbe = -5
  integer, parameter :: chf_pbe0 = -6
  integer, parameter :: chf_lcwpbe = -7
  integer, parameter :: chf_pw86 = -8
  integer, parameter :: chf_b971 = -9

  ! overloaded functions
  interface realloc
     module procedure realloc1r
     module procedure realloc2r
     module procedure realloc1i
  end interface

contains
  subroutine error (routine,message,errortype)
    character*(*), intent(in) :: routine !< routine calling the error
    character*(*), intent(in) :: message !< the message
    integer, intent(in) :: errortype !< fatal, warning or info
    
    character*(20)      chtype

    if (errortype == 2) then
       chtype='ERROR'
    else if (errortype == 1) then
       chtype='WARNING'
    else if (errortype == 0) then
       chtype='COMMENT'
    else
       chtype='UNKNOWN'
    endif
    write (iout,100) trim(adjustl(chtype)),trim(adjustl(routine)),&
       trim(adjustl(message))

    if (errortype == 2) stop 1

100 format (A,"(",A,"): ",A)
    
  end subroutine error

  function elem2z(s) result(z)
    character*2, intent(in) :: s
    integer :: z

    integer :: i
    
    z = -1
    do i = 0, 103
       if (s == ptable(i)) then
          z = i
          return
       endif
    enddo

  end function elem2z

  function z2nr(z)
    integer, intent(in) :: z
    integer :: z2nr

    z2nr = 40
    if (z > 2) z2nr = 60
    if (z > 10) z2nr = 80
    if (z > 18) z2nr = 100
    if (z > 36) z2nr = 120
    if (z > 54) z2nr = 140
    if (z > 86) z2nr = 160

  endfunction z2nr

  function z2nang(z)
    integer, intent(in) :: z
    integer :: z2nang

    z2nang = 194

  endfunction z2nang

  function frevol(z)

    integer, intent(in) :: z
    real*8 :: frevol

    ! DKH LSDA/UGBS Free Atomic Volumes
    real*8, parameter :: frevol0(0:103) = (/&
       0.d0,&
       9.194D0,   4.481D0,  91.957D0,  61.357D0,  49.813D0,  36.728D0,&
       27.633D0,  23.517D0,  19.322D0,  15.950D0, 109.359D0, 103.064D0,&
       120.419D0, 104.229D0,  86.782D0,  77.133D0,  66.372D0,  57.336D0,&
       203.093D0, 212.202D0, 183.101D0, 162.278D0, 143.250D0, 108.209D0,&
       123.098D0, 105.735D0,  92.944D0,  83.794D0,  75.750D0,  81.177D0,&
       118.371D0, 116.334D0, 107.474D0, 103.221D0,  95.111D0,  87.605D0,&
       248.772D0, 273.748D0, 249.211D0, 223.801D0, 175.809D0, 156.831D0,&
       160.042D0, 136.654D0, 127.754D0,  97.024D0, 112.778D0, 121.627D0,&
       167.906D0, 172.030D0, 165.500D0, 163.038D0, 153.972D0, 146.069D0,&
       341.992D0, 385.767D0, 343.377D0, 350.338D0, 334.905D0, 322.164D0,&
       310.337D0, 299.537D0, 289.567D0, 216.147D0, 268.910D0, 259.838D0,&
       251.293D0, 243.174D0, 235.453D0, 228.284D0, 229.617D0, 209.971D0,&
       197.541D0, 183.236D0, 174.685D0, 164.139D0, 150.441D0, 135.765D0,&
       125.297D0, 131.258D0, 185.769D0, 195.671D0, 193.036D0, 189.142D0,&
       185.919D0, 181.089D0, 357.787D0, 407.283D0, 383.053D0, 362.099D0,&
       346.565D0, 332.462D0, 319.591D0, 308.095D0, 297.358D0, 300.572D0,&
       275.792D0, 266.317D0, 257.429D0, 209.687D0, 203.250D0, 230.248D0,&
       236.878D0/)

    real*8, parameter :: frevol_blyp(0:36) = (/0.0d0,&
       8.6751280810827840d0,  4.3645522950717863d0,  89.719664495180297d0,  61.278566735200307d0,&
       49.519428604126382d0,  36.855287686102294d0,  27.995151669578650d0,  23.764998306645893d0,&
       19.585085705364346d0,  16.198873185715303d0,  110.72206109235110d0,  104.99418885629647d0,&
       124.86364864987824d0,  105.67572783021383d0,  87.717325899499158d0,  77.848560257666222d0,&
       67.289916054772661d0,  58.134006480572936d0,                0.00d0,                0.00d0,&
       189.14940377667708d0,  171.31844286164150d0,  149.09546424400813d0,  114.14072711487664d0,&
       128.27850058556785d0,  108.89501751641639d0,  100.22570628848071d0,  88.525764030777225d0,&
       79.682125997948589d0,  86.320056678306727d0,  125.54483633211143d0,  121.00196086145128d0,&
       110.76811978260324d0,  106.26389443012715d0,  98.394240990660435d0,  90.352795412830417d0/)

    real*8, parameter :: frevol_b3lyp(0:36) = (/0.0d0,&
       8.3489695345990036d0,   4.2304165146306838d0,   88.786443805705602d0,   60.659572325858221d0,&
       48.180704821245364d0,   35.739786376347816d0,   27.113474668623077d0,   22.954305454378311d0,&
       18.919035719008878d0,   15.664743495672592d0,   110.10781957063872d0,   104.58688460317792d0,&
       121.82866689851529d0,   103.64265930616017d0,   86.233229483159292d0,   76.534831202674539d0,&
       66.225085507692555d0,   57.258468294377920d0,                 0.00d0,                 0.00d0,&
       191.24458275864666d0,   171.33960226014415d0,   154.31690117090523d0,   111.84266101395626d0,&
       128.53304731251799d0,   116.41853334783917d0,   108.02725265724374d0,   84.734127226732753d0,&
       78.821628181721010d0,   86.690015512660722d0,   122.56838049698425d0,   118.57117842661357d0,&
       108.88205342457658d0,   104.55842043337798d0,   96.968864870029904d0,   89.143550757884739d0/)

    real*8, parameter :: frevol_bhahlyp(0:36) = (/0.0d0,&
       8.0162383356457330d0,   4.0855044443453972d0,   89.521965961222222d0,   60.766620056468859d0,&
       47.163104244372654d0,   34.732519196898075d0,   26.261248810982419d0,   22.163686460916928d0,&
       18.257837007647932d0,   15.123544024501175d0,   112.24236986192338d0,   105.88927909724653d0,&
       120.53914101629064d0,   102.42697418515390d0,   85.163324240312747d0,   75.599161377353852d0,&
       65.445419322876759d0,   56.592207744291457d0,                 0.00d0,                 0.00d0,&
       196.29478033696057d0,   175.14878735114101d0,   157.82044436715444d0,   112.66706636630363d0,&
       131.50432127575891d0,   120.17186617063507d0,   110.69809115434559d0,   102.41457915920718d0,&
       79.806952905523545d0,   88.942968328671427d0,   121.67064092786063d0,   117.25525753962938d0,&
       107.64004292512693d0,   103.49115496709287d0,   96.074020468035101d0,   88.360967365230920d0/)

    real*8, parameter :: frevol_camb3lyp(0:36) = (/0.0d0,&
       8.4912418798953802d0,   4.3046720809866850d0,   88.495476866973974d0,   60.708687086264462d0,&
       48.017293545760708d0,   35.742021952558758d0,   27.186861831258817d0,   23.059619196158092d0,&
       19.022193061849432d0,   15.756088275099161d0,   109.34523733398444d0,   104.23708524052515d0,&
       119.90584483487682d0,   102.80245055238821d0,   85.960831224028368d0,   76.454433439613283d0,&
       66.279667949988067d0,   57.375114690853394d0,                 0.00d0,                 0.00d0,&
       192.25186364675062d0,   171.60008159670065d0,   154.54149488373386d0,   111.51614466860299d0,&
       105.82575187101907d0,   97.244815362081795d0,   90.518519096216096d0,   84.233565981799387d0,&
       78.777041862773274d0,   86.830269580599790d0,   119.93150277142236d0,   117.12848860140892d0,&
       108.22197759832967d0,   104.14881532532733d0,   96.808752555226278d0,   89.137842527944699d0/)

    real*8, parameter :: frevol_pbe(0:36) = (/0.0d0,&
       8.7017290470668396d0,   4.3452296013273557d0,   90.345016391875745d0,   60.583880555189161d0,&
       49.405610809309636d0,   36.699238724745918d0,   27.804307905805462d0,   23.504820577931341d0,&
       19.366722293791256d0,   16.015015785319534d0,   114.11568983081936d0,   105.34255389919765d0,&
       122.27570960010875d0,   104.11673088180522d0,   86.659855994802655d0,   76.696371875869517d0,&
       66.330422698625398d0,   57.338280708600728d0,                 0.00d0,                 0.00d0,&
       187.56578096056526d0,   171.16842233126320d0,   148.57152896051596d0,   113.40726679106993d0,&
       128.41370162442536d0,   107.47276564091520d0,   98.987585702398064d0,   87.019526141682050d0,&
       79.950858621796073d0,   86.179295605944361d0,   123.32536942981970d0,   119.49220381812802d0,&
       109.53573383239170d0,   104.78641436128312d0,   97.036678428323626d0,   89.131992369076343d0/)

    real*8, parameter :: frevol_pbe0(0:36) = (/0.0d0,&
       8.2794385587230224d0,   4.1765568461439084d0,   89.990651247904850d0,   60.150518670639258d0,&
       47.925007649610208d0,   35.403450375407488d0,   26.774856262986901d0,   22.577665436425793d0,&
       18.604506038051770d0,   15.403636840460491d0,   114.72060124074002d0,   105.56153583642701d0,&
       119.94555610703479d0,   102.21839478488732d0,   85.124898059747593d0,   75.344227406670839d0,&
       65.219744182377752d0,   56.417325659383977d0,                 0.00d0,                 0.00d0,&
       191.49634784581932d0,   172.37802114577727d0,   155.38093772078227d0,   111.48097347813784d0,&
       129.94413257650865d0,   118.24195915775856d0,   108.60919379569880d0,   84.698903894685841d0,&
       79.451081251238904d0,   87.077874282230383d0,   120.92043418719176d0,   117.15417864570686d0,&
       107.60403177731271d0,   103.07063942096924d0,   95.595468315759462d0,   87.905274745875047d0/)

    real*8, parameter :: frevol_lcwpbe(0:36) = (/0.0d0,&
       8.2370522934694321d0,   4.3223022069392556d0,   88.889190621676747d0,   59.167955275706255d0,&
       46.644645536860530d0,   35.000149688018325d0,   26.874237675680991d0,   22.830136179756057d0,&
       18.961662156609091d0,   15.787003768893198d0,   115.52305481049332d0,   105.20996979820804d0,&
       115.58130215151486d0,   99.382061772188109d0,   83.591029347688959d0,   74.254876662748089d0,&
       64.635098656590586d0,   56.191387396031978d0,                 0.00d0,                 0.00d0,&
       191.95897600535326d0,   172.74995562085155d0,   155.27567955463721d0,   110.25556775344671d0,&
       129.26382250173134d0,   117.09487172302605d0,   107.19200454429361d0,   83.698218802709803d0,&
       78.438819467992630d0,   85.426026058281309d0,   114.73764666891768d0,   113.36022260917842d0,&
       105.50701995208520d0,   101.36644049619710d0,   94.506386358808882d0,   87.304323578968550d0/)

    real*8, parameter :: frevol_pw86(0:36) = (/0.0d0,&
       8.5848597505149957d0,   4.3045044427345758d0,   89.105017427126995d0,   60.116125959883448d0,&
       48.948725444378958d0,   36.639492704666679d0,   27.890608359059900d0,   23.482178111604007d0,&
       19.392669723193279d0,   16.068196774925887d0,   110.77812603458771d0,   103.32122266742829d0,&
       121.98526328640968d0,   104.32311632196750d0,   86.866106676280893d0,   76.813580132172405d0,&
       66.479441471513695d0,   57.487390010861731d0,                 0.00d0,                 0.00d0,&
       185.73638891654375d0,   168.40972266496894d0,   146.57709921981228d0,   114.30230293458976d0,&
       127.03123778878130d0,   106.92921677417986d0,   79.779742298798681d0,   87.455617557953332d0,&
       79.471676692123495d0,   85.082884646536542d0,   123.88868121486288d0,   120.63195902241362d0,&
       110.68925500617975d0,   105.58567728669152d0,   97.798158973601787d0,   89.844087877827207d0/)

    real*8, parameter :: frevol_b971(0:36) = (/0.0d0,&
       8.2753044447676203d0,   4.1538172892410463d0,   92.548449604186061d0,   60.079247546006165d0,&
       47.862796246981773d0,   35.305970222823603d0,   26.699111937155749d0,   22.644191719940501d0,&
       18.675265336738804d0,   15.464799224948079d0,   121.31651007378332d0,   105.68066738843501d0,&
       120.26430745377735d0,   101.96967180988391d0,   85.008007953828127d0,   75.622028007253363d0,&
       65.525564639136547d0,   56.692633782240122d0,                 0.00d0,                 0.00d0,&
       193.61923913260267d0,   171.80907721991764d0,   154.23278648930369d0,   114.13224357718498d0,&
       107.38843940948567d0,   98.695645988131020d0,   92.056435539203278d0,   86.394941236711290d0,&
       80.059364904746786d0,   87.457830726555457d0,   121.47718849519879d0,   116.93338080131248d0,&
       107.07435287066555d0,   103.20705595068554d0,   95.879065901986593d0,   88.146215884084469d0/)

    real*8 :: rchf

    ! gaussian basis set, %HF=(0,25,50,75,100) for periods 1 and 2.
    real*8 frevol1(5,0:10)
    save frevol1
    data frevol1/&
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
       8.7017290470668396d0, 8.2794385587230224d0, 7.9275669093485428d0, 7.6272753143131924d0, 7.3655545253236543d0,&
       4.3452296013273557d0, 4.1765568461439084d0, 4.0288815322803657d0, 3.8977496795733644d0, 3.7798775917054876d0,&
       90.345016391875745d0, 89.990651247904850d0, 89.843447659370099d0, 89.862545631289919d0, 90.018126381136568d0,&
       60.583880555189161d0, 60.150518670639258d0, 59.789689626639294d0, 59.490664095467629d0, 59.244839599587813d0,&
       49.405610809309636d0, 47.925007649610208d0, 46.767078359393196d0, 45.827336436528093d0, 45.040918655916336d0,&
       36.699238724745918d0, 35.403450375407488d0, 34.355914497332748d0, 33.486552145008638d0, 32.747999094847692d0,&
       27.804307905805462d0, 26.774856262986901d0, 25.931460984272665d0, 25.224788496044660d0, 24.620417659921344d0,&
       23.504820577931341d0, 22.577665436425793d0, 21.811051051552955d0, 21.165993022839672d0, 20.614016545011395d0,&
       19.366722293791256d0, 18.604506038051770d0, 17.967854713848929d0, 17.427825507512161d0, 16.962774503019283d0,&
       16.015015785319534d0, 15.403636840460491d0, 14.886731288556328d0, 14.443878154969715d0, 14.059464430879629d0/

    ! pure GGA
    if (abs(chf) < 1d-10) then
       if (z > 10) then
          frevol = frevol0(z)
       else
          frevol = frevol1(1,z)
       endif
       return
    endif

    ! special cases
    if (chf < 0d0 .and. (z<19.or.z>20.and.z<37)) then ! up to Kr except K and Ca
       select case(nint(chf))
       case(chf_blyp) 
          frevol = frevol_blyp(z)
       case(chf_b3lyp) 
          frevol = frevol_b3lyp(z)
       case(chf_bhahlyp) 
          frevol = frevol_bhahlyp(z)
       case(chf_camb3lyp) 
          frevol = frevol_camb3lyp(z)
       case(chf_pbe) 
          frevol = frevol_pbe(z)
       case(chf_pbe0) 
          frevol = frevol_pbe0(z)
       case(chf_lcwpbe) 
          frevol = frevol_lcwpbe(z)
       case(chf_pw86) 
          frevol = frevol_pw86(z)
       case(chf_b971) 
          frevol = frevol_b971(z)
       case default
          call error("frevol","unknown functional",2)
       end select
    else
       ! general hybrid
       if (chf < 0d0) then
          select case(nint(chf))
          case(chf_blyp) 
             rchf = 0d0
          case(chf_b3lyp) 
             rchf = 0.2d0
          case(chf_bhahlyp) 
             rchf = 0.5d0
          case(chf_camb3lyp) 
             rchf = 0.2d0
          case(chf_pbe) 
             rchf = 0.0d0
          case(chf_pbe0) 
             rchf = 0.25d0
          case(chf_lcwpbe) 
             rchf = 0.25d0
          case(chf_pw86) 
             rchf = 0.0d0
          case(chf_b971) 
             rchf = 0.21d0
          case default
             call error("frevol","unknown functional",2)
          end select
       else
          rchf = chf
       endif

       if (z > 10) then
          frevol = frevol0(z)
       else
          if (rchf < 0.25d0) then
             frevol = frevol1(1,z) + (frevol1(2,z)-frevol1(1,z)) / 0.25d0 * (rchf-0d0)
          elseif (rchf < 0.50d0) then
             frevol = frevol1(2,z) + (frevol1(3,z)-frevol1(2,z)) / 0.25d0 * (rchf-0.25d0)
          elseif (rchf < 0.75d0) then
             frevol = frevol1(3,z) + (frevol1(4,z)-frevol1(3,z)) / 0.25d0 * (rchf-0.50d0)
          else
             frevol = frevol1(4,z) + (frevol1(5,z)-frevol1(4,z)) / 0.25d0 * (rchf-0.75d0)
          endif
       endif
    endif

  endfunction frevol

 !> Convert string to lowercase except where quoted
  !> string and lower will return the same
  function lower (string)
    
    character*(*), intent(in) :: string !< Input string, same as the function result on output
    character*(len(string)) :: lower

    integer :: i, iadd

    lower = string
    iadd = ichar('A') - ichar('a')
    i = 1
    do while (i<=len(string))
       if (string(i:i).ge.'A' .and. string(i:i).le.'Z') then
          lower(i:i) = char (ichar(string(i:i)) - iadd)
       endif
       i = i + 1
    enddo
  end function lower

  !> Get integer value from input text. If a valid integer is not
  !> found, then return .false.
  logical function isinteger (ival,line,lp)
    
    character*(mline), intent(in) :: line !< Input string
    integer, intent(inout) :: lp !< Pointer to current position on string
    integer, intent(out) :: ival !< Integer value read, same as function result

    integer :: i

    isinteger = .false.
    do while (line(lp:lp).eq.blank)
       lp=lp+1
       if (lp > mline) return
    enddo
    i=lp
    if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i=i+1
    if (isdigit(line(i:i))) then 
       do while (isdigit(line(i:i)))
          i=i+1
       enddo
       if (line(i:i) .eq. blank .or. line(i:i) .eq. null .or. line(i:i).eq.newline ) then
          ival=atoi(line(lp:i-1))
          lp = i
          isinteger=.true.
       else
          ival=0
          isinteger=.false.
       endif
    else
       ival=0
       isinteger=.false.
    endif
  end function isinteger

  !> Get a real number from line and sets rval to it.
  !> If a valid real number is not found, isreal returns .false.
  logical function isreal (rval, line, lp)
    
    character*(*), intent(in) :: line !< Input string
    integer, intent(inout) :: lp !< Pointer to current position on string
    real*8, intent(out) :: rval !< Real value read

    character*(mline) :: dumchar
    integer           :: tp, i
    character*(1)     :: ch
    logical           :: matched

    isreal = .false.
    do while (line(lp:lp).eq.blank)
       lp = lp + 1
       if (lp > mline) return
    end do

    i = lp
    if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1 
    if (isdigit(line(i:i))) then 
       do while (isdigit(line(i:i)))
          i = i + 1
       enddo
       if (line(i:i) .eq. '.') then 
          i = i + 1
          do while (isdigit(line(i:i)))
             i = i + 1
          enddo
       endif
       matched = .true.
    else if (line(i:i) .eq. '.') then 
       i = i + 1
       if (isdigit(line(i:i))) then 
          do while (isdigit(line(i:i)))
             i = i + 1
          enddo
          matched = .true.
       else
          matched = .false.
       endif
    else
       matched = .false.
    endif

    !.....get optional exponent
    tp = i - 1
    if (matched) then 
       if (line(i:i)=='e' .or. line(i:i)=='E' .or. line(i:i)=='d' .or. line(i:i)=='D'.or.&
           line(i:i)=='-' .or. line(i:i)=='+') then 
          i = i + 1
          if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1 
          if (isdigit (line(i:i))) then 
             do while (isdigit(line(i:i)))
                i = i + 1
             enddo
             if (index(blank//','//null//newline, line(i:i)) .gt. 0) then
                dumchar=line(lp:i-1)//null
                rval = atof (dumchar)
                lp = i
             else
                matched = .false.
                rval = 0d0
             endif
          else 
             matched = .false.
          endif
       else
          if (index(blank//','//null//newline, line(i:i)) .gt. 0) then
             dumchar=line(lp:tp)//null
             rval = atof(dumchar)
             lp = i
          else
             matched = .false.
             rval = 0d0
          endif
       endif
    else
       rval = 0d0
    endif
    !
    isreal = matched
    return

  end function isreal

  !> Convert ascii string to integer (private).
  integer function atoi (string)
    
    character*(*), intent(in) :: string !< Input string

    integer           i, sign

    atoi=0
    i=1
    do while (string(i:i).eq.blank)
       i=i+1
    end do
    sign=1
    if(string(i:i) .eq. '+' .or. string(i:i) .eq. '-') then 
       if (string(i:i) .eq. '-') sign=-1
       i=i+1
    endif
    do while (isdigit(string(i:i)))
       atoi=10*atoi+ichar(string(i:i))-ichar('0')
       i=i+1
       if (i > len(string)) then
          exit
       end if
    end do
    atoi=atoi*sign
    return

  end function atoi

  !> Return true if c is a digit (private).
  logical function isdigit (c)
    
    character*(1), intent(in) :: c !< Is the character c a digit?

    isdigit = c.ge.'0' .and. c.le.'9'

  end function isdigit

  !> Convert ascii string to real (private).
  real*8 function atof (str)

    character*(*),intent(in) :: str !< Input string

    real*8, parameter :: ten=10d0

    real*8 val,power
    integer exponent,sign,esign,i

    sign=1
    val=0d0
    power=1d0
    exponent=0
    esign=1
    i=1
    do while (str(i:i) .eq. blank)
       i=i+1 
    end do
    if (str (i:i) .eq. '+' .or. str(i:i) .eq.'-') then 
       if (str(i:i) .eq. '-') sign=-1
       i=i+1
    endif
    do while (isdigit(str(i:i)))
       val=ten*val+ichar(str(i:i))-ichar('0')
       i=i+1
    enddo
    if (str(i:i) .eq. '.') then 
       i=i+1
       do while (isdigit(str(i:i)))
          val=ten*val+ichar(str(i:i))-ichar('0')
          i=i+1
          power=power*ten
       enddo
    endif
    if (str(i:i).eq.'e' .or. str(i:i).eq.'E' .or. str(i:i).eq.'d' .or. str(i:i).eq.'D') then
       i=i+1
       if (str (i:i) .eq. '+' .or. str(i:i) .eq.'-') then 
          if (str(i:i) .eq. '-') esign=-1
          i=i+1
       endif
       do while (isdigit(str(i:i)))
          exponent=10*exponent+ichar(str(i:i))-ichar('0')
          i=i+1
       end do
    elseif (str(i:i).eq.'-' .or. str(i:i).eq.'+') then
       esign = 1
       if (str(i:i) .eq. '-') esign=-1
       i = i + 1 
       do while (isdigit(str(i:i)))
          exponent=10*exponent+ichar(str(i:i))-ichar('0')
          i=i+1
       end do
    endif

    atof=(sign*val/power)*ten**(esign*exponent)
    return
  end function atof

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc1r(a,nnew)

    real*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    real*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call error('realloc1r','array not allocated',2)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1r

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc2r(a,n1,n2)

    real*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    real*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call error('realloc2r','array not allocated',2)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    allocate(temp(n1,n2))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2r

  !> Adapt the size of an allocatable 1D integer array
  subroutine realloc1i(a,nnew)

    integer, intent(inout), allocatable :: a(:) !< Input array, integer, 1D
    integer, intent(in) :: nnew !< New dimension
    
    integer, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call error('realloc1i','array not allocated',2)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1i

end module param
