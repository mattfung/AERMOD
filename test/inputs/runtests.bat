REM Standard Test Cases for AERMOD model, 
REM dated 14134, May 14, 2014

echo %date% %time% > start.txt

copy aertest.inp  aermod.inp
aermod
copy aermod.out   aertest.out

copy aertest_evt.inp  aermod.inp
aermod
copy aermod.out       aertest_evt.out

copy allsrcs.inp  aermod.inp
aermod
copy aermod.out   allsrcs.out

copy allsrcs_evt.inp aermod.inp
aermod
copy aermod.out      allsrcs_evt.out

copy capped.inp  aermod.inp
aermod
copy aermod.out  capped.out

copy capped_evt.inp aermod.inp
aermod
copy aermod.out     capped_evt.out

copy capped_nostd.inp aermod.inp
aermod
copy aermod.out       capped_nostd.out

copy capped_nostd_evt.inp aermod.inp
aermod
copy aermod.out           capped_nostd_evt.out

copy testgas.inp  aermod.inp
aermod
copy aermod.out   testgas.out

copy testgas2.inp aermod.inp
aermod
copy aermod.out   testgas2.out

copy testgas2_evt.inp aermod.inp
aermod
copy aermod.out       testgas2_evt.out

copy testpart.inp aermod.inp
aermod
copy aermod.out   testpart.out

copy testprt2.inp aermod.inp
aermod
copy aermod.out   testprt2.out

copy testprt2_evt.inp aermod.inp
aermod
copy aermod.out       testprt2_evt.out

copy openpits.inp   aermod.inp
aermod
copy aermod.out   openpits.out

copy openpits_evt.inp aermod.inp
aermod
copy aermod.out     openpits_evt.out

copy lovett.inp   aermod.inp
aermod
copy aermod.out   lovett.out

copy lovett_evt.inp aermod.inp
aermod
copy aermod.out     lovett_evt.out

copy flatelev.inp   aermod.inp
aermod
copy aermod.out   flatelev.out

copy flatelev_evt.inp aermod.inp
aermod
copy aermod.out     flatelev_evt.out

copy mcr.inp  aermod.inp
aermod
copy aermod.out   mcr.out

copy mcr_evt.inp  aermod.inp
aermod
copy aermod.out   mcr_evt.out

copy multurb.inp  aermod.inp
aermod
copy aermod.out   multurb.out

copy multurb_evt.inp aermod.inp
aermod
copy aermod.out      multurb_evt.out

copy hrdow.inp    aermod.inp
aermod
copy aermod.out   hrdow.out

copy hrdow_evt.inp aermod.inp
aermod
copy aermod.out    hrdow_evt.out

copy surfcoal.inp aermod.inp
aermod
copy aermod.out   surfcoal.out

copy surfcoal_evt.inp aermod.inp
aermod
copy aermod.out       surfcoal_evt.out

copy olm.inp      aermod.inp
aermod
copy aermod.out   olm.out

copy olm_evt.inp  aermod.inp
aermod
copy aermod.out   olm_evt.out

copy olmgrp.inp   aermod.inp
aermod
copy aermod.out   olmgrp.out

copy olmgrp_evt.inp  aermod.inp
aermod
copy aermod.out      olmgrp_evt.out

copy pvmrm.inp    aermod.inp
aermod
copy aermod.out   pvmrm.out

copy pvmrm_evt.inp aermod.inp
aermod
copy aermod.out    pvmrm_evt.out

copy psdcred.inp  aermod.inp
aermod
copy aermod.out   psdcred.out

copy scimtest.inp  aermod.inp
aermod
copy aermod.out    scimtest.out


rem MULTYEAR PM10 test case

copy testpm10_1986.inp aermod.inp
aermod
copy aermod.out testpm10_1986.out

copy testpm10_1987.inp aermod.inp
aermod
copy aermod.out testpm10_1987.out

copy testpm10_1988.inp aermod.inp
aermod
copy aermod.out testpm10_1988.out

copy testpm10_1989.inp aermod.inp
aermod
copy aermod.out testpm10_1989.out

copy testpm10_1990.inp aermod.inp
aermod
copy aermod.out testpm10_1990.out

rem PM10 EVT test with MULTYEAR option

copy testpm10_multyr_evt.inp aermod.inp
aermod
copy aermod.out testpm10_multyr_evt.out


rem PM10 test case with 5-yr data file

copy testpm10.inp aermod.inp
aermod
copy aermod.out testpm10.out

rem PM10 EVT test with 5-yr option

copy testpm10_evt.inp aermod.inp
aermod
copy aermod.out testpm10_evt.out


rem PM2.5 test case with 5-yr data file

copy testpm25.inp aermod.inp
aermod
copy aermod.out   testpm25.out

rem PM2.5 EVT test case with 5-yr data file

copy testpm25_evt.inp aermod.inp
aermod
copy aermod.out       testpm25_evt.out

echo %date% %time% > end.txt