print('------ FUNCTION SSP.MODEL ------')
print('--------------------------------')
print('PURPOSE:')
print('   Create SSP spectra')
print('CALLING SEQUENCE:')
print('   ssp.model(feh, afe, age, imf, slope, fwhm, dl, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe,')
print('             n_ms, n_rg, parfile, lmin, lmax)')
print('INPUTS:')
print('   feh   = iron abundance [Fe/H]')
print('   afe   = [alpha/Fe]')
print('   age   = age of the isochrone (Gyr)')
print('   imf   = IMF (Kroupa, Salpeter or Unimodal)')
print('   slope = IMF slope (if IMF = Unimodal)')
print('   fwhm  = spectral resolution (default = 0.2 A)')
print('   dl    = sampling delta lambda (default = 0.1 A/pixel)')
print('   CFe   = [C/Fe] (set to 0.0 if omitted, i.e., default [C/Fe] = solar)')
print('   NFe   = [N/Fe] (set to 0.0 if omitted, i.e., default [N/Fe] = solar)')
print('   OFe   = [O/Fe] (set to [alpha/Fe] if omitted, i.e., default [O/Fe] = afe)')
print('   MgFe  = [Mg/Fe] (set to [alpha/Fe] if omitted, i.e., default [Mg/Fe] = afe)')
print('   SiFe  = [Si/Fe] (set to [alpha/Fe] if omitted, i.e., default [Si/Fe] = afe)')
print('   CaFe  = [Ca/Fe] (set to [alpha/Fe] if omitted, i.e., default [Ca/Fe] = afe)')
print('   TiFe  = [Ti/Fe] (set to a[alpha/Fe]fe if omitted, i.e., default [Ti/Fe] = afe)')
print('   NaFe  = [Na/Fe] (set to 0.0 if omitted, i.e., default [Na/Fe] = solar)')
print('   n_ms  = number of main sequence stars')
print('   n_rg  = number of red giant stars')
print('   parfile = file with list of stellar parameters (used only if n_ms/n_rg are not specified)')
print('   lmin  = lower lambda')
print('   lmax  = upper lambda')
print('OUTPUT:')
print('   SSP spectra in folder ./SSP_Spectra/' )
print('   logfile: SSP_model.log' )
print('REQUIRED FUNCTIONS:')
print('   nulbadgrade (fortran code)')
print('   set.stpars.filename (R function, defined in stpars.R)')
print('   get.tracks          (R function, defined in stpars.R)')
print('   set.stspec.filename (R function, defined in pfant12.R)')
print('EXAMPLES:')
print("   ssp.model(0.2, 0.2, 8, 'Kroupa', n_ms = 9, n_rg = 6)")
print("   ssp.model(0.2, 0.2, 8, 'Unimodal', slope = 2.3, n_ms = 9, n_rg = 6)")


####################################################################
# FUNCTION SSP.MODEL
####################################################################
ssp.model <- function(feh, afe, age, imf, slope, fwhm, dl, CFe, NFe,
                      OFe, MgFe, SiFe, CaFe, TiFe, NaFe, n_ms, n_rg, parfile, lmin, lmax){

  #---------------------------------
  # CHECKING INPUTS
  #---------------------------------
  if(missing(feh)){print('Need to specify [Fe/H] (feh=)')}
  if(missing(afe)){print('Need to specify [alpha/Fe] (afe=)')}
  if(missing(age)){print('Need to specify population age in Gyr (age=)')}
  if(missing(imf)){print('Need to specify IMF (Salpeter, Kroupa, Unimodal)')}
  imf = tolower(imf)
  if(imf == 'unimodal' & missing(slope)){
    print('IMF = unimodal, you need to specify slope')
  }
  
  if(missing(n_ms) & missing(n_rg) & missing(parfile)){
    print('FILE WITH STELLAR PARAMETERS --> ???')
    print(' Need to specify n_ms & n_rg OR parfile')
  }else{
    if(missing(n_ms) == F & missing(n_rg) == F){
      pars.file = set.stpars.filename(n_ms, n_rg, feh, afe, age)      
    }else{
      pars.file = parfile
    }
  }
  
  #---------------------------------
  # SETTING DEFAULT VALUES
  #---------------------------------
  if(missing(CFe)){CFe = 0.0}
  if(missing(NFe)){NFe = 0.0}
  if(missing(OFe)){OFe = afe}
  if(missing(MgFe)){MgFe = afe}
  if(missing(SiFe)){SiFe = afe}
  if(missing(CaFe)){CaFe = afe}
  if(missing(TiFe)){TiFe = afe}
  if(missing(NaFe)){NaFe = 0.0}
  if(missing(fwhm)){fwhm = 0.2}
  if(missing(dl)){dl = 0.1}

  #---------------------------------
  # OUTPUT FILE NAMES
  #---------------------------------
  file_ssp <- set.ssp.filename(feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe, age, imf, slope)
  file_fig = paste(file_ssp, '.jpg', sep = '')
  logfile = paste(file_ssp, '.log', sep = '')  

  #---------------------------------
  # WRITE INFO IN LOGFILE
  #---------------------------------
  write(paste('----------------------- INPUT PARAMETERS -----------------------'), file = logfile, append = F)
  write(paste('----------------------------------------------------------------'), file = logfile, append = T)
  
  write(paste('System time: ', Sys.time()), file = logfile, append = T)
  
  temp <- sprintf('%-26s%8.2f%4s', 'Isochrone age: ', age, ' Gyr')
  write(temp, file = logfile, append = T)
  temp <- sprintf('%-26s%8.2f%2s', 'FWHM: ', fwhm, ' A')
  write(temp, file = logfile, append = T)
  if(imf == 'unimodal'){
    temp <- sprintf('%-26s%8s%8s%8.2f', 'IMF: ', imf, ' slope: ', slope)
  }else{
    temp <- sprintf('%-26s%8s', 'IMF: ', imf)
  }
  write(temp, file = logfile, append = T)
  
  write(paste('Abundances:'), file = logfile, append = T)
  temp <- sprintf('%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s', '[Fe/H]', '[a/Fe]', 
                  '[C/Fe]', '[N/Fe]', '[O/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Ti/Fe]', '[Na/Fe]')
  write(temp, file = logfile, append = T)
  
  temp <- sprintf('%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f', 
                  feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
  write(temp, file = logfile, append = T)
  
  #---------------------------------
  # CONSTANTS
  #---------------------------------
  Lsun = 3.839e33  # erg s-1
  
  G = 6.67300e-11  # m3 kg-1 s-2
  Msun = 1.9891e33  # g
  sb = 5.67e-5 # erg cm-2 s-1 K-4
  Rsun = sqrt(Lsun / (4 * pi * sb * 5777**4)) #  = 6.955e10 cm
  aSun = 4 * pi * Rsun**2  # cm2
  
  #---------------------------------
  # READ INPUTS FILE
  #---------------------------------
  t = read.table(pars.file)
  nstars = nrow(t)
  Teffs = t$V1
  loggs = t$V2
  masses = t$V3 
  lumis = 10**t$V4 # Lsun
  radii = sqrt(lumis * Lsun / (4 * pi * sb * Teffs**4)) # cm
  area = 4 * pi * (radii)**2      # cm2

  write(paste('Reading stellar parameters from: ', pars.file), file = logfile, append = T)
  write(paste('----------------------------------------------------------------\n'), file = logfile, append = T)

  
  #---------------------------------
  # MASS BINS
  #---------------------------------
  logL = log10(lumis)
  logL_breaks = logL
  nn = length(lumis)
  
  lum_breaks = lumis
  
  for(i in 1:(nn + 1)){
    if(i == 1){lum_breaks[i] = lumis[1] }
    if(i > 1 & i < (nn + 1)){lum_breaks[i] = mean(c(lumis[i - 1], lumis[i]))}
    if(i == (nn + 1)){lum_breaks[i] = lumis[nn]}
  }
 
  iso.file <- get.tracks(feh, afe, age)
  t = read.table(iso.file)
  iso.mass = t$V2
  iso.logL = t$V5
  iso.L = 10**iso.logL
  mass_breaks <- interp(iso.logL, log10(lum_breaks), iso.mass) 
  mass_breaks[1] = 0.08
  mass_breaks[(nn+1)] = masses[nn]
 
  # Correct mass_breaks if interpolation fails for stars close to the turn-off
  for(i in 1:nn){
    if(mass_breaks[i] > mass_breaks[i + 1] | mass_breaks[i] > masses[i]){
      mass_breaks[i] = masses[i] - (masses[i] - masses[i-1])/2
      logL_breaks[i] = logL[i] - (logL[i] - logL[i-1])/2
    }
  }
  
  #---------------------------------
  # SSP SPECTRA
  #---------------------------------
  fraction_M = masses; fraction_M[] = 0
  fraction_L = masses; fraction_L[] = 0
  sum_flux = masses; sum_flux[] = 0
  L_corr =  masses; L_corr[] = 0

  write(paste('SSP divided in: ', nstars, 'mass bins'), file = logfile, ncolumns = 100, append = T)  

  write('Information about each bin:', file = logfile, ncolumns = 100, append = T)
  write('----------------------', file = logfile, ncolumns = 1, append = T)
  data <- sprintf('%5s%10s%10s%10s%10s%7s%7s%10s%10s%10s%10s%10s%10s%10s', 
                  'i', 'logL_i', 'logL_f', 'Mass_i', 'Mass_f',
                  'Teff', 'logg', 'Mass', 'logLstar', 'Radius', 'logLbin', 'Lcorr', 'frac_M', 'frac_L')
  write(data, file = logfile, ncolumns = 100, append = T)
  write('    -----------------------------------------------------------------------------------------------------------------------------',
        file = logfile, ncolumns = 100, append = T)
  for(i in 1:nstars){
    file_flux <- set.stspec.filename(feh, afe, lmin, lmax, Teffs[i], 
                                     loggs[i], CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
    
    #---------------------------------
    # CONVOLUTION STELLAR SPECTRA
    #---------------------------------
    system('rm -f st_spectra')
    system(paste('cp ', file_flux, ' temp.in'))
    data = c('./nulbadgrade << DONE', 'temp.in', 'F', 'T', 'st_spectra', dl, 'T', paste(fwhm), 'DONE')
    write(data, file = 'nulbadgrade.sh', sep = '\n', append = F)
    system('chmod 777 nulbadgrade.sh')
    system('./nulbadgrade.sh')
    
    #---------------------------------
    # READ CONVOLVED STELLAR SPECTRA
    #---------------------------------
    tt = read.table('st_spectra')
    lambda = tt$V1
    units = 1e-8 * 4 * pi   # erg / (s cm2 cm ster) --> erg / (s cm2 A)
    # Scale flux according to the star surface area and normalize by Lsun
    st_flux = tt$V2 * units * area[i] / Lsun  # erg / (s cm2 A) --> Lsun / A
    
    if(i == 1){
      ssp_flux = st_flux; ssp_flux[] = 0
    }

    # Number of stars formed between mass_breaks[i] and mass_breaks[i+1] (per unit of mass):
    #  sum(phi_m * dm) = phi_mdm
    dm_break = mass_breaks[i+1] - mass_breaks[i]
    dm = dm_break/1000
    if(dm > 0.01){dm = 0.001}
    masses_temp = seq(mass_breaks[i], mass_breaks[i+1], dm)
    iso.L_temp <- 10**interp(iso.mass, masses_temp, iso.logL) 
    phi_mdm = 0
    mass_bin = 0
    L_bin = 0
    L_star_bin = 0
    
    #---------------------------------
    # IMF: KROUPA
    #---------------------------------
    if(imf == 'kroupa'){
      for(j in 1:length(masses_temp)){
        cc = 1
        if(masses_temp[j] >= 1.0  & masses_temp[j] <  100){x = 2.7}
        if(masses_temp[j] >= 0.5  & masses_temp[j] <    1){
          x = 2.3
          cc = ((1.0**(-2.7)) / 1.0**(-2.3))
        }
        if(masses_temp[j] >= 0.08 & masses_temp[j] <  0.5){
          x = 1.3
          cc = ((0.5**(-2.3)) / 0.5**(-1.3)) * ((1.0**(-2.7)) / 1.0**(-2.3))
        }
        if(masses_temp[j] >= 0.01 & masses_temp[j] < 0.08){
          x = 0.3
          cc = (0.08**(-1.3)) / 0.08**(-0.3) * 
            ((0.5**(-2.3)) / 0.5**(-1.3)) * ((1.0**(-2.7)) / 1.0**(-2.3))
        }
        phi_m = masses_temp[j]**(-x) * cc
        phi_mdm = phi_mdm +  phi_m * dm
        mass_bin = mass_bin + masses_temp[j] * phi_m * dm
        L_bin = L_bin + iso.L_temp[j] * phi_m * dm
        L_star_bin = L_star_bin + lumis[i] * phi_m * dm
      }
    }
    
    #---------------------------------
    # IMF: SALPETER
    #---------------------------------
    if(imf == 'salpeter'){
      x = 2.3
      for(j in 1:length(masses_temp)){
        phi_m = masses_temp[j]**(-x)
        phi_mdm = phi_mdm +  phi_m * dm
        mass_bin = mass_bin + masses_temp[j] * phi_m * dm
        L_bin = L_bin + iso.L_temp[j] * phi_m * dm
        L_star_bin = L_star_bin + lumis[i] * phi_m * dm
      }
    }
    
    #---------------------------------
    # IMF: UNIMODAL (SLOPE=x)
    #---------------------------------
    if(imf == 'unimodal'){
      x = slope
      for(j in 1:length(masses_temp)){
        phi_m = masses_temp[j]**(-x)
        phi_mdm = phi_mdm + phi_m * dm
        mass_bin = mass_bin + masses_temp[j] * phi_m * dm
        L_bin = L_bin + iso.L_temp[j] * phi_m * dm
        L_star_bin = L_star_bin + lumis[i] * phi_m * dm
      }  
    }

    #---------------------------------
    # SUM(S * PHI * DM)
    #---------------------------------
    L_corr[i] = L_bin / L_star_bin
    ssp_flux = ssp_flux + st_flux * L_corr[i] * phi_mdm
    fraction_M[i] = mass_bin
    sum_flux[i] = sum(st_flux) * (max(lambda) -  min(lambda))
    fraction_L[i] = L_bin
  }
  
  # normalization: sum(mass * phi_m * dm) = 1 Msun
  L_SSP = sum(fraction_L) 
  fraction_L = fraction_L / L_SSP
  ssp_flux = ssp_flux / sum(fraction_M)   
  fraction_M = fraction_M / sum(fraction_M)

  #---------------------------------
  # WRITE INFO IN LOGFILE
  #---------------------------------

  for(i in 1:nstars){
    data <- sprintf('%5i%10.4f%10.4f%10.4f%10.4f%7.0f%7.2f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.2f', 
                     i, log10(lum_breaks[i]), log10(lum_breaks[i+1]), mass_breaks[i], mass_breaks[i+1], Teffs[i], 
                     loggs[i], masses[i], log10(lumis[i]), radii[i] / Rsun, log10(L_SSP * fraction_L[i]),
                     L_corr[i], fraction_M[i]*100, fraction_L[i]*100)
    write(data, file = logfile, ncolumns = 100, append = T)
  }

  write('----------------------', file = logfile, ncolumns = 1, append = T)
  data = c('logL_i, logL_f   -> lower and upper luminosity limits of the bin (log[L/Lsun])', 
           'Mass_i, Mass_f   -> lower and upper mass limits of the bin (Msun)',
           'Teff, logg       -> stelar parameter of the synthetic stellar spectrum representing the stars within bin',
           'Mass, logLstar   -> mass and luminosity of the star with Teff, logg',
           'Radius           -> radius of a star with Teff, logg (Rsun, Radius = sqrt[L / (4 * pi * sb * Teffs**4)])',
           'logLbin          -> total luminosity of the bin (log[L/Lsun])',
           'Lcorr            -> luminosity correction factor taking into account the variation of stellar luminosities',
           '                      within bin ( = sum[L(m)*phi_m*dm] / sum[Lstar*phi_m*dm])',
           'frac_M           -> fraction of mass in bin i (%, frac_M = (Mbin/M_SSP)*100)',
           'frac_L           -> fraction of light in bin i (%, frac_L = (Lbin/L_SSP)*100)')

  write(data, file = logfile, ncolumns = 1, append = T)
   
  data <- sprintf('%7.2f%13.5e', lambda[i], ssp_flux[i])
  write(data, file = 'SSP_spectra', ncolumns = 2)
  for(i in 2:length(lambda)){
    data <- sprintf('%7.2f%13.5e', lambda[i], ssp_flux[i])
    write(data, file = 'SSP_spectra', append = T, ncolumns = 100)
  }

  system(paste('cp SSP_spectra ', file_ssp))
  
  write(paste('\n    SSP spectra -->  ', file_ssp), file = logfile, ncolumns = 100, append = T)
  write(paste('----------------------------------------------------------------\n'), file = logfile, append = T)
  
  #---------------------------------
  # PLOT SSP SPECTRA
  #---------------------------------
  t = read.table(file_ssp)
  jpeg(file_fig, width = 800, height = 500)
  plot(t$V1, t$V2, type = 'l', xlab = expression(paste(lambda, ' (', ring(A), ')')), ylab = 'Flux')
  dev.off()
}


##################################
# FUNCTION INTERP
##################################
interp <- function(xin, xout, yin){
  Ninterpol = length(xout)
  yout <- vector(mode = 'numeric', length = Ninterpol)
  
  for(k in 1:Ninterpol){
    t = xin[xin < xout[k]]
    tind = length(t)
    
    if(tind <= 1){tind = 2}
    if(tind >= length(xin)){tind = length(xin) - 1}
    t1 = xin[tind - 1]
    t2 = xin[tind]
    t3 = xin[tind + 1]
    tx = xout[k]
    
    A = (tx - t1) / (t3 - t1)
    B = (tx - t2) / (t3 - t2)
    C = (tx - t3) / (t2 - t1)
    D = (tx - t1) / (t3 - t2)
    E = (tx - t2) / (t3 - t1)
    
    yout[k] = yin[tind+1] * A * B - yin[tind] * D * C + yin[tind-1] * E * C
  } 
  return(yout)
}


##################################
# FUNCTION SET.SSP.FILENAME
##################################
set.ssp.filename <- function(feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe, age, imf, slope){
  #---------------------------------
  # SETTING DEFAULT VALUES
  #---------------------------------
  if(missing(CFe)){CFe = 0.0}
  if(missing(NFe)){NFe = 0.0}
  if(missing(OFe)){OFe = afe}
  if(missing(MgFe)){MgFe = afe}
  if(missing(SiFe)){SiFe = afe}
  if(missing(CaFe)){CaFe = afe}
  if(missing(TiFe)){TiFe = afe}
  if(missing(NaFe)){NaFe = 0.0}
  imf = tolower(imf)
  
  temp = c(feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
  temp_s = vector()
  for(j in 1:length(temp)){
    if(temp[j] < 0){
      temp_s[j] = sprintf('%1s%4.2f', '-', abs(temp[j]))
    }else{
      temp_s[j] = sprintf('%1s%4.2f', '+', temp[j])
    }
  }
  
  if(age < 10){age_s = sprintf('%1s%3.1f', '0', age)}else{age_s = sprintf('%4.1f', age)}
  
  if(imf == 'unimodal'){
    slope_s = sprintf('%4.2f', slope)
    file_ssp = paste('./SSP_Spectra/SSP_Fe', temp_s[1],  '_a', temp_s[2], '_C', temp_s[3], '_N', temp_s[4],
                      '_O', temp_s[5],'_Mg', temp_s[6], '_Si', temp_s[7],'_Ca', temp_s[8],'_Ti', temp_s[9], 
                     '_Na', temp_s[10], '_age', age_s, '_slope', slope_s, sep = '')
  }else{
    file_ssp = paste('./SSP_Spectra/SSP_Fe', temp_s[1],  '_a', temp_s[2], '_C', temp_s[3], '_N', temp_s[4],
                      '_O', temp_s[5],'_Mg', temp_s[6], '_Si', temp_s[7],'_Ca', temp_s[8],'_Ti', temp_s[9], 
                     '_Na', temp_s[10], '_age', age_s, '_', imf, sep = '')
  }
}







