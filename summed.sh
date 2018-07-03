name="modeles"
rm -f ${name}.mod
rm -f ${name}mod.log
touch ${name}.mod
touch ${name}mod.log
metl=" -0.5 "
Teff=" 3267 "
logg=" 0.03 "
echo '----------- STARTING THE MEGA-CYCLE -----------'
  lin=0
  lin=`expr $lin + 1`
  ./innewmarcs2 << DONE
  '${name}.mod' '${name}.dat'
  T
  ${Teff} ${logg} ${metl} newmarcs
  ${lin}
  F
  DONE
