print('------ FUNCTION STPARS ------')
print('-----------------------------')
print('PURPOSE:')
print('   Get stellar parameters from stellar evolutionary tracks')
print('CALLING SEQUENCE:')
print('   stpars(n_ms, n_rg, feh, afe, age, fig)')
print('INPUTS:')
print('   n_ms = number of desired main sequence stars')
print('   n_rg = number of desired red giant stars')
print('   feh  = iron abundance [Fe/H] (Available [Fe/H] range is -0.5 -> 0.5')
print('   afe  = [alpha/Fe] (available values -0.2, 0.0, +0.2, +0.4, +0.6, +0.8')
print('                    [alpha/Fe] = +0.4, +0.6, +0.8 -> available only for [Fe/H] <=0)')
print('   age  = age of the population (in Gyr; available ages: see manual)')
print('   fig  = logical; set to T to plot the isochrone')
print('OUTPUT:')
print('   text file containing the stellar parameters (Teff, logg, Mass, logL)' )
print('     of (n_ms + n_rg + 1) stars' )
print('REQUIRED SCRIPTS:')
print('   ./atm_models/iso_interp_feh (fortran code)')
print('   ./atm_models/isolf_split (fortran code)')
print('EXAMPLE:')
print('   stpars(9, 6, 0.2, 0.2, 8)')

####################################################################
# FUNCTION STPARS
####################################################################
stpars <- function(n_ms, n_rg, feh, afe, age, fig){
  #---------------------------------
  # SET OUTPUT FILENAME
  #---------------------------------
  output.file <- set.stpars.filename(n_ms, n_rg, feh, afe, age)
  
  #---------------------------------
  # ISOCHRONE [FE/H]-INTERPOLATION
  #---------------------------------
  iso.file <- get.tracks(feh, afe, age)
  
  #---------------------------------
  # READ ISOCHRONE
  #---------------------------------
  t = read.table(iso.file)
  iso.teff = 10**t$V3
  iso.logg = t$V4
  iso.mass = t$V2
  iso.logL = t$V5

  #---------------------------------
  # GET STELLAR PARAMETERS
  #---------------------------------
  logg_lim = max(iso.logg[iso.teff == max(iso.teff)])
  ms = iso.logg >= logg_lim
  rg = iso.logg < logg_lim
  min_teff_ms = min(iso.teff[ms])
  min_teff_rg = min(iso.teff[rg])
  
  delta_teff_ms = (max(iso.teff) - min_teff_ms) / n_ms
  delta_teff_rg = (max(iso.teff) - min_teff_rg) / n_rg
  
  teff_grid = c(1:(n_ms + n_rg + 1))
  
  for(i in 1:n_ms){
    teff_grid[i] = min_teff_ms + (i - 1) * delta_teff_ms
  }
  
  teff_grid[n_ms + 1] = max(iso.teff)
  j = 0
  for(i in (n_ms + n_rg + 1):(n_ms + 2)){
    teff_grid[i] = min_teff_rg + j * delta_teff_rg
    j = j + 1
  }
  
  index_lim = which(teff_grid == max(iso.teff))
  logg_grid = teff_grid
  mass_grid = teff_grid
  lumi_grid = teff_grid
  
  for(i in 1: length(teff_grid)){
    if(i <= index_lim){
      xx.xx = iso.logg >= logg_lim
    }else{
      xx.xx = iso.logg < logg_lim
    }
    temp = abs(iso.teff[xx.xx] - teff_grid[i] )
    teff_grid[i] = iso.teff[xx.xx][temp == min(temp)][1]
    logg_grid[i] = iso.logg[xx.xx][temp == min(temp)][1]
    mass_grid[i] = iso.mass[xx.xx][temp == min(temp)][1]
    lumi_grid[i] = iso.logL[xx.xx][temp == min(temp)][1]
  }

  iso.teff.grid = teff_grid
  iso.logg.grid = logg_grid
  iso.mass.grid = mass_grid
  iso.logL.grid = lumi_grid
  
  #---------------------------------
  # WRITE OUTPUT FILE
  #---------------------------------
  data <- sprintf('%10s%10s%12s%10s', '# Teff/K', 'logg', 'Mass/Msun',  'logL/Lsun')
  write(data, file = paste(output.file, sep =''), ncolumns = 100)
  for(i in 1:(n_ms + n_rg + 1)){
    data <- sprintf('%10.0f%10.2f%12.7f%10.4f', 
                    iso.teff.grid[i], iso.logg.grid[i], iso.mass.grid[i], iso.logL.grid[i])
    write(data, file = paste(output.file, sep =''), 
          ncolumns = 100, append = T)
  }
  
  #---------------------------------
  # PLOT ISOCHRONE
  #---------------------------------
  if(missing(fig) == F & fig != F){
    par(mar = c(5, 6, 4, 3))
    plot(1, 1, xlab = 'Temperature (K)', ylab = 'log g (dex)', pch = 19, cex = 1.5, 
         ylim = c(5.5, -0.2), xlim = c(8000, 2000), cex.lab = 1.6, cex.axis = 1.6)
    lines(iso.teff, iso.logg, lwd = 2, col = 'red')
    points(iso.teff.grid, iso.logg.grid, col = 'red', pch = 19, cex = 1.4)
    legend(8000, 0, c(paste('Isochrone [Fe/H]=', feh, ', [a/Fe]=', afe, 
                            ', Age=', age, ' Gyr', sep = ''), 
                      'Selected stellar parameters'), bty = 'n', cex = 0.8,
           lty = c(1, 1), col = c('red', 'red'), pch = c(19, 19), pt.cex = c(0, 1), 
           lwd = c(2, 0))
  }
  
  
  #---------------------------------
  # PRINT SOME INFORMATION
  #---------------------------------
  print(paste('Isochrone: [Fe/H] = ', feh))
  print(paste('           [a/Fe] = ', afe))
  print(paste('           Age    = ', age, 'Gyr')) 
  print(paste('STELLAR PARAMETERS OF: ', n_ms + n_rg + 1, 'STARS IN THE'))
  print(paste('OUTPUT FILE: ', output.file))
  
}

##################################
# FUNCTION GET.TRACKS
##################################
get.tracks <- function(feh, afe, age){
  afes = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8)
  temp = sort(abs(afes - afe), index.return = T)
  
  script.file = 'iso.sh'
  write('./iso_interp_feh << DONE', file = script.file, append = F)
  write('11', file = script.file, append = T)
  write('1', file = script.file, append = T)
  write(temp$ix[1], file = script.file, append = T)
  write(feh, file = script.file, append = T)
  write('isotemp', file = script.file, append = T)
  write('DONE', file = script.file, append = T)
  
  write('#', file = script.file, append = T)
  
  write('./isolf_split << DONE', file = script.file, append = T)
  write('isotemp', file = script.file, append = T)
  write('DONE', file = script.file, append = T)
  
  system('chmod 777 iso.sh')
  system('mv iso.sh ./DARTMOUTH/iso.sh')
  setwd('./DARTMOUTH/')
  system('./iso.sh > temp.txt')
  setwd('../')
  
  if(age < 10){
    iso.file.out = paste('./DARTMOUTH/a0', age * 1000, 'isotemp', sep = '')
  }else{
    iso.file.out = paste('./DARTMOUTH/a', age * 1000, 'isotemp', sep = '')
  }
  return(iso.file.out)
}

##################################
# FUNCTION SET.STPARS.FILENAME
##################################
set.stpars.filename <- function(n_ms, n_rg, feh, afe, age){
  if(feh < 0){feh_s = sprintf('%1s%4.2f', '-', abs(feh))}else{feh_s = sprintf('%1s%4.2f', '+', feh)}
  if(afe < 0){afe_s = sprintf('%1s%4.2f', '-', abs(afe))}else{afe_s = sprintf('%1s%4.2f', '+', afe)}
  if(age < 10){age_s = sprintf('%1s%3.1f', '0', age)}else{age_s = sprintf('%4.1f', age)}
  if(n_ms < 10){n_ms_s = sprintf('%1s%1i', '0', n_ms)}else{n_ms_s = sprintf('%2i', n_ms)}
  if(n_rg < 10){n_rg_s = sprintf('%1s%1i', '0', n_rg)}else{n_rg_s = sprintf('%2i', n_rg)}
  file = paste('./Stellar_pars/stpars_', 'fe', feh_s, 'a', afe_s, 'age', age_s, 
                      'ms', n_ms_s, 'rg', n_rg_s, '.dat', sep = '')
  
  return(file)
}


