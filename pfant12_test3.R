print('------ FUNCTION PFANT12 ------')
print('------------------------------')
print('PURPOSE:')
print('   Calculate synthetic stellar spectra using the PFANT code')
print('CALLING SEQUENCE:')
print('   pfant12(parfile, feh, afe, lmin, lmax, age, vt, fwhm, dl, CFe, NFe,')
print('           OFe, MgFe, SiFe, CaFe, TiFe, NaFe, n_ms, n_rg, parfile)')
print('INPUTS:')
print('   feh  = iron abundance [Fe/H]')
print('   afe  = [alpha/Fe]')
print('   lmin = lower lambda')
print('   lmax = upper lambda')
print('   age  = age of the isochrone (Gyr)')
print('   vt   = microturbulence velocity (default = 2.0 km/s)')
print('   fwhm = spectral resolution (default = 0.2 A)')
print('   dl   = sampling delta lambda (default = 0.1 A/pixel)')
print('   CFe  = [C/Fe] (set to 0.0 if omitted, i.e., default [C/Fe] = solar)')
print('   NFe  = [N/Fe] (set to 0.0 if omitted, i.e., default [N/Fe] = solar)')
print('   OFe  = [O/Fe] (set to [alpha/Fe] if omitted, i.e., default [O/Fe] = afe)')
print('   MgFe = [Mg/Fe] (set to [alpha/Fe] if omitted, i.e., default [Mg/Fe] = afe)')
print('   SiFe = [Si/Fe] (set to [alpha/Fe] if omitted, i.e., default [Si/Fe] = afe)')
print('   CaFe = [Ca/Fe] (set to [alpha/Fe] if omitted, i.e., default [Ca/Fe] = afe)')
print('   TiFe = [Ti/Fe] (set to a[alpha/Fe]fe if omitted, i.e., default [Ti/Fe] = afe)')
print('   NaFe = [Na/Fe] (set to 0.0 if omitted, i.e., default [Na/Fe] = solar)')
print('   n_ms = number of main sequence stars')
print('   n_rg = number of red giant stars')
print('   parfile = file with list of stellar parameters (used only if n_ms/n_rg are not specified)')
print('OUTPUT:')
print('   Synthetic stellar spectra in folder ./Stellar_Spectra/' )
print('   logfile: pfant12.log' )
print('REQUIRED SCRIPTS:')
print('   pfantgrade (fortran code)')
print('   nulbadegrade (fortran code)')
print('   hydro2 (fortran code)')
print('   innewmarcs2 (fortran code)')
print('   set.stpars.filename (R function, defined in stpars.R)')
print('EXAMPLE:')
print('   pfant12(0.2, 0.2, 6000, 6500, 8, n_ms = 9, n_rg = 6)')

# Defining PATH to PYFANT
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/home/tatiana/PFANT/fortran/bin/", sep=":"))

####################################################################
# FUNCTION PFANT12
####################################################################
pfant12 <- function(feh, afe, lmin, lmax, age, vt, fwhm, dl, CFe, NFe,
                    OFe, MgFe, SiFe, CaFe, TiFe, NaFe, n_ms, n_rg, parfile){
  logfile = 'pfant12.log'
  #---------------------------------
  # CHECKING INPUTS
  #---------------------------------
  if(missing(feh)){print('Need to specify [Fe/H] (feh=)')}
  if(missing(afe)){print('Need to specify [alpha/Fe] (afe=)')}
  if(missing(lmin) | missing(lmax)){print('Need to specify wavelength range (lmin=, lmax=)')}
  if(missing(age)){print('Need to specify population age in Gyr (age=)')}
  
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
  if(missing(vt)){vt = 2.0}
  
  
  #---------------------------------
  # WRITE INFO IN LOGFILE
  #---------------------------------
  write(paste('----------------------- INPUT PARAMETERS -----------------------'), file = logfile, append = T)
  write(paste('----------------------------------------------------------------'), file = logfile, append = T)
  
  write(paste('System time: ', Sys.time()), file = logfile, append = T)
  
  temp <- sprintf('%-26s%8.2f%4s', 'Isochrone age: ', age, ' Gyr')
  write(temp, file = logfile, append = T)
  temp <- sprintf('%-26s%8.0f%3s%6.0f%2s', 'Wavelength limits: ', lmin, ' - ', lmax, ' A')
  write(temp, file = logfile, append = T)
  temp <- sprintf('%-26s%8.2f%5s', 'Microturbulence velocity: ', vt, ' km/s')
  write(temp, file = logfile, append = T)
  temp <- sprintf('%-26s%8.2f%2s', 'FWHM: ', fwhm, ' A')
  write(temp, file = logfile, append = T)
  
  write(paste('Abundances:'), file = logfile, append = T)
  temp <- sprintf('%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s', '[Fe/H]', '[a/Fe]', 
                  '[C/Fe]', '[N/Fe]', '[O/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Ti/Fe]', '[Na/Fe]')
  write(temp, file = logfile, append = T)
  
  temp <- sprintf('%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f', 
                  feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
  write(temp, file = logfile, append = T)
  
  #---------------------------------
  # READ INPUT FILE
  #---------------------------------
  t = read.table(pars.file)
  Teffs = t$V1
  loggs = t$V2
  nstars = nrow(t)
  
  write(paste('Reading stellar parameters from: ', pars.file), file = logfile, append = T)
  write(paste('----------------------------------------------------------------\n'), file = logfile, append = T)
  write(paste('----------------------------------------------------------------\n'), file = logfile, append = T)
  
  #---------------------------------
  # CALCULATE SYNTHETIC SPECTRA
  #---------------------------------
  for(i in c(1:nstars)){
    #---------------------------------
    # SETTING FILE NAME
    #---------------------------------
    file_flux <- set.stspec.filename(feh, afe, lmin, lmax, Teffs[i], 
                                     loggs[i], CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
    
    file_fig  = paste(file_flux, '.jpg', sep = '') 
    
    file_flux_cont = paste(file_flux, '_cont', sep = '')
    file_flux_norm = paste(file_flux, '_norm', sep = '')
    file_flux_convol = paste(file_flux, '_convol', sep = '')
    
    #---------------------------------
    # CHECK IF SPECTRUM EXISTS
    #---------------------------------
    flag.missing = 0
    flag.missing = file.access(file_flux, mode = 0)
    if(flag.missing == 0){
      write(paste('Synthetic stellar spectra: (Teff, logg) = (', Teffs[i], ', ', loggs[i], ')', 
                  ' already exists --> skipping', sep = ''), file = logfile, append = T)
      write(paste('     file: ', file_flux, sep = ''), file = logfile, append = T)
      write(paste('----------------------------------------------------------------\n'), file = logfile, append = T)
    }else{
      write(paste('Calculating synthetic stellar spectra: (Teff, logg) = (', Teffs[i], ', ', loggs[i], ')', sep = ''),
            file = logfile, append = T)
      write(paste('Starting at: ', Sys.time()), file = logfile, append = T)
      ptm <- proc.time()
      
      #write.abonds.dissoc(CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
      
      vt_s = vt
      Teff_s = Teffs[i]
      logg_s = loggs[i]
      met_grade_s = feh
      
      #---------------------------------
      # INTERPOLATION MARCS MODELS
      #---------------------------------
      system('innewmarcs --allow T --flprefix session-12/flux --fn_abonds session-12/abonds.dat --fn_dissoc session-12/dissoc.dat --fn_hmap session-12/hmap.dat --fn_main session-12/main.dat --fn_modeles session-12/modeles.mod --fn_opa session-12/modeles.opa --fn_progress session-12/progress.txt  --opa F')
      
      #---------------------------------
      # CALCULATE HYDROGEN LINES
      #---------------------------------
      system('hydro2 --allow T --flprefix session-12/flux --fn_abonds session-12/abonds.dat --fn_dissoc session-12/dissoc.dat --fn_hmap session-12/hmap.dat --fn_main session-12/main.dat --fn_modeles session-12/modeles.mod --fn_opa session-12/modeles.opa --fn_progress session-12/progress.txt --opa F')
      
      #---------------------------------
      # RUN THE MAIN CODE (PFANT)
      #---------------------------------
      system('pfant --allow T --flprefix session-12/flux --fn_abonds session-12/abonds.dat --fn_dissoc session-12/dissoc.dat --fn_hmap session-12/hmap.dat --fn_main session-12/main.dat --fn_modeles session-12/modeles.mod --fn_opa session-12/modeles.opa --fn_progress session-12/progress.txt --opa F')
      
      #---------------------------------
      # CONVOLUTION
      #---------------------------------
      system('nulbad --allow T --flprefix session-12/flux --fn_abonds session-12/abonds.dat --fn_dissoc session-12/dissoc.dat --fn_hmap session-12/hmap.dat --fn_main session-12/main.dat --fn_modeles session-12/modeles.mod --fn_opa session-12/modeles.opa --fn_progress session-12/progress.txt --opa F')
      
      #---------------------------------
      # WRITE STELLAR SPECTRA
      #---------------------------------
      system(paste('cp session-12/flux.spec ', file_flux, sep = '')) 
      system(paste('cp session-12/flux.cont ', file_flux_cont, sep = '')) 
      system(paste('cp session-12/flux.norm ', file_flux_norm, sep = '')) 
      system(paste('cp session-12/flux.norm.nulbad.',fwhm,'00 ',file_flux_convol, sep = '')) 
      
      #---------------------------------
      # PLOT STELLAR SPECTRA
      #---------------------------------
      t = read.table(file_flux_convol)
      jpeg(file_fig, width = 800, height = 500)
      plot(t$V1, t$V2, type = 'l', xlab = expression(paste(lambda, ' (', ring(A), ')')), 
           ylab = 'Flux')
      dev.off()
      
      #---------------------------------
      # WRITE INFO IN LOGFILE
      #---------------------------------
      time <- proc.time() - ptm 
      write(paste('Finishing at: ', Sys.time()), file = logfile, append = T)
      write(sprintf('%14s%10.0f%4s%-10.3f%5s', 'Elapsed time: ', time[3], ' s (', time[3]/60, ' min)'), 
            file = logfile, append = T)
      write(paste('----------------------------------------------------------------\n'), file = logfile, append = T)
    }
  }  
}


