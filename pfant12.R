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
        
        write.abonds.dissoc(CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
        
        vt_s = vt
        Teff_s = Teffs[i]
        logg_s = loggs[i]
        met_grade_s = feh
        
        #---------------------------------
        # INTERPOLATION MARCS MODELS
        #---------------------------------
        create.model(Teff_s, logg_s, met_grade_s, logfile)
        
        #---------------------------------
        # CALCULATE HYDROGEN LINES
        #---------------------------------
        system('cp absoru2_hydro2.dat absoru2.dat')
        data = c('./hydro2 << DONE', '1', 'F', '2 3 6562.817 10.199 2442.326', 'T', vt_s, 
                 Teff_s, logg_s, met_grade_s, '1', 'thalpha', 'F', 'DONE', '######',
                 './hydro2 << DONE', '1', 'F', '2 4 4861.342 10.199  249.628', 'T', vt_s, 
                 Teff_s, logg_s, met_grade_s, '1', 'thbeta', 'F', 'DONE', '######',
                 './hydro2 << DONE', '1', 'F', '2 5 4340.475 10.199  74.4776', 'T', vt_s, 
                 Teff_s, logg_s, met_grade_s, '1', 'thgamma', 'F', 'DONE', '######',
                 './hydro2 << DONE', '1', 'F', '2 6 4101.742 10.199  32.8903', 'T', vt_s, 
                 Teff_s, logg_s, met_grade_s, '1', 'thdelta', 'F', 'DONE', '######',
                 './hydro2 << DONE', '1', 'F', '2 7 3970.070 10.199  17.72642', 'T', vt_s, 
                 Teff_s, logg_s, met_grade_s, '1', 'thepsilon', 'F', 'DONE')
        write(data, file = 'hydro2.sh', sep = '\n', append = F)
        system('chmod 777 hydro2.sh')
        system('./hydro2.sh')
        system('cp absoru2_pfant.dat absoru2.dat')
        system('cp thepsilon th1')
        
        #---------------------------------
        # WRITE MAIN.DAT
        #---------------------------------
        l1 = 'syn'
        l2 = paste('T ', dl, ' 5.0 1.   .12')
        l3 = vt_s
        l4 = paste(Teff_s, logg_s, met_grade_s, '0.1 1')
        l5 = 'F  1.'
        l6 = met_grade_s
        l7 = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
        l8 = 'flux'
        l9 = paste(lmin, lmax, ' 50')
        l10 = 'th1\nth1\nth1\nth1\nth1\nthepsilon\nthdelta\nthgamma\nthbeta\nthalpha'
        
        data = c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10)
        write(data, file = 'main.dat', sep = '\n', append = F)
        
        #---------------------------------
        # RUN THE MAIN CODE (PFANT)
        #---------------------------------
        system('./pfantgrade')
        
        #---------------------------------
        # CONVOLUTION
        #---------------------------------
        system('rm -f int')
        data = c('./nulbadgrade << DONE', 'spec.flux', 'F', 'T', 'int', dl, 'T', paste(fwhm), 'DONE')
        write(data, file = 'nulbadgrade.sh', sep = '\n', append = F)
        system('chmod 777 nulbadgrade.sh')
        system('./nulbadgrade.sh')
        
        #---------------------------------
        # WRITE STELLAR SPECTRA
        #---------------------------------
        system(paste('cp spec.flux ', file_flux, sep = '')) 
        system(paste('cp cont.flux ', file_flux_cont, sep = '')) 
        system(paste('cp norm.flux ', file_flux_norm, sep = '')) 
        system(paste('cp int ', file_flux_convol, sep = '')) 
        
        
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

##################################
# FUNCTION WRITE.ABONDS.DISSOC
##################################
write.abonds.dissoc <- function(CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe){
  # DISSOC DATA
  #---------------------------------
  data1.dissoc = "   18  100    0.005   -1.0"
  data3.dissoc = c("AN      1.28051E01 -8.27934E00 6.41622E-02 -7.36267E-03 3.46663E-04 2 91 71",
                   "CN      1.28051E01 -8.27934E00 6.41622E-02 -7.36267E-03 3.46663E-04 2 61 71",
                   "CAH     1.13401E01 -3.01442E00 4.23487E-01 -6.14674E-02 3.16392E-03 2201 11",
                   "MGO     1.17018E01 -5.03261E00 2.96408E-01 -4.28111E-02 2.20232E-03 2121 81",
                   "TIO     1.33981E01 -8.59562E00 4.08726E-01 -5.79369E-02 2.92873E-03 2221 81",
                   "MGH     1.12853E01 -2.71637E00 1.96585E-01 -2.73103E-02 1.38164E-03 2121 11",
                   "AC      1.28038E01 -6.51780E00 9.77186E-02 -1.27393E-02 6.26035E-04 2 91 61",
                   "AA      1.28038E01 -6.51780E00 9.77186E-02 -1.27393E-02 6.26035E-04 1 92   ",
                   "CC      1.28038E01 -6.51780E00 9.77186E-02 -1.27393E-02 6.26035E-04 1 62   ",
                   "HOH     2.54204E01 -1.05223E01 1.69394E-01 -1.83684E-02 8.17296E-04 2 12 81",
                   "ALH     1.21913E01 -3.76361E00 2.55568E-01 -3.72612E-02 1.94061E-03 2131 11",
                   "ALO     1.27393E01 -5.25336E00 1.82177E-01 -2.57927E-02 1.31850E-03 2131 81",
                   "AO      1.38200E01 -1.17953E01 1.72167E-01 -2.28885E-02 1.13491E-03 2 91 81",
                   "CO      1.38200E01 -1.17953E01 1.72167E-01 -2.28885E-02 1.13491E-03 2 61 81",
                   "AH      1.21355E01 -4.07599E00 1.27676E-01 -1.54727E-02 7.26615E-04 2 91 11",
                   "CH      1.21355E01 -4.07599E00 1.27676E-01 -1.54727E-02 7.26615E-04 2 61 11",
                   "AP      1.24988E01 -6.47155E00 1.23609E-01 -1.74113E-02 8.97797E-04 2 91151",
                   "CP      1.24988E01 -6.47155E00 1.23609E-01 -1.74113E-02 8.97797E-04 2 61151",
                   "AS      1.34357E01 -8.55736E00 1.87544E-01 -2.55069E-02 1.27346E-03 2 91161",
                   "CS      1.34357E01 -8.55736E00 1.87544E-01 -2.55069E-02 1.27346E-03 2 61161",
                   "NH      1.20327E01 -3.84349E00 1.36286E-01 -1.66427E-02 7.86913E-04 2 71 11",
                   "NO      1.28310E01 -7.19644E00 1.73495E-01 -2.30652E-02 1.13799E-03 2 71 81",
                   "NN      1.35903E01 -1.05855E01 2.20671E-01 -2.99975E-02 1.49927E-03 1 72   ",
                   "OH      1.23710E01 -5.05783E00 1.38217E-01 -1.65474E-02 7.72245E-04 2 81 11",
                   "OO      1.32282E01 -5.51807E00 6.99354E-02 -8.15109E-03 3.79699E-04 1 82   ",
                   "NAH     1.14155E01 -2.75079E00 1.96460E-01 -2.73828E-02 1.38445E-03 2111 11",
                   "SIF     1.24272E01 -5.83726E00 1.66854E-01 -2.27883E-02 1.13711E-03 2141 91",
                   "SIH     1.18522E01 -3.74185E00 1.59988E-01 -2.06292E-02 9.98967E-04 2141 11",
                   "SIO     1.34132E01 -8.87098E00 1.50424E-01 -1.95811E-02 9.48283E-04 2141 81",
                   "SIN     1.23989E01 -5.48756E00 9.53008E-02 -1.33693E-02 6.93956E-04 2141 71",
                   "PH      1.20802E01 -4.64388E00 3.41137E-01 -4.87959E-02 2.50567E-03 2151 11",
                   "CAO     1.22598E01 -6.05249E00 5.82844E-01 -8.58050E-02 4.44251E-03 2201 81",
                   "SCO     1.37467E01 -8.64196E00 4.80722E-01 -6.96697E-02 3.57468E-03 2211 81",
                   "FEO     1.29845E01 -5.33182E00 3.17452E-01 -4.45649E-02 2.25240E-03 2261 81",
                   "NIH     1.25203E01 -3.43307E00 1.96303E-01 -2.53774E-02 1.24117E-03 2281 11",
                   "CUH     1.21518E01 -3.91900E00 3.09765E-01 -4.33886E-02 2.19329E-03 2291 11",
                   "CUO     1.21756E01 -4.28270E00 2.04849E-01 -2.82166E-02 1.40747E-03 2291 81",
                   "HH      1.27388E01 -5.11717E00 1.25720E-01 -1.41494E-02 6.30214E-04 1 12   ", '')
  data3.dissoc <- sprintf("%75s", data3.dissoc)
    
  # ABUNDANCES: [X/H] -> STELLAR IRON ABUNDANCE [Fe/H] ALREADY ADDED
  #    (VALUES IN DISSOC AND ABONDS FILES ARE (X/H) - [Fe/H])
  # INITIAL ABUNDANCES (SOLAR)
  abonds.ini.file = 'abonds_Sun.dat'
  dissoc.ini.file = 'dissoc_Sun.dat'
  
  # READ INITIAL ABUNDANCES
  #---------------------------------
  # ABONDS.DAT
  t = read.table(abonds.ini.file, fill = T, as.is = 1)
  NABOND = length(t$V1) - 2
  
  elems = t$V1[1:NABOND]
  elems <- sprintf("%2s", elems)
  
  abonds.ini = t$V2[1:NABOND]
  abonds.fin = abonds.ini
  for(k in 1:NABOND){
    if(elems[k] == ' C'){abonds.fin[k] = CFe + abonds.ini[k]}
    if(elems[k] == ' N'){abonds.fin[k] = NFe + abonds.ini[k]}
    if(elems[k] == ' O'){abonds.fin[k] = OFe + abonds.ini[k]}
    if(elems[k] == 'MG'){abonds.fin[k] = MgFe + abonds.ini[k]}
    if(elems[k] == 'SI'){abonds.fin[k] = SiFe + abonds.ini[k]}
    if(elems[k] == 'CA'){abonds.fin[k] = CaFe + abonds.ini[k]}
    if(elems[k] == 'TI'){abonds.fin[k] = TiFe + abonds.ini[k]}
    if(elems[k] == 'NA'){abonds.fin[k] = NaFe + abonds.ini[k]}
  } # END for(k in 1:NABOND)

  # DISSOC.DAT
  t = read.table(dissoc.ini.file, fill = T, as.is = 1)
  elems.dissoc  = t$V1
  elems.dissoc <- sprintf("%2s", elems.dissoc)
  dissoc.c2     = t$V2
  dissoc.c3     = t$V3
  dissoc.c4     = t$V4
  dissoc.c5     = t$V5
  abonds.dissoc = t$V6
  
  NABOND.DISSOC = length(elems.dissoc)
  abonds.dissoc.fin = abonds.dissoc
 
  # DISSOC.DAT = ABONDS.DAT - 12
  for(k1 in 1:NABOND){
    for(k2 in 1:NABOND.DISSOC){
      if(elems.dissoc[k2] == elems[k1]){
        abonds.dissoc.fin[k2] = abonds.ini[k1]-12
      }
    }
  }

  #---------------------------------
  # WRITE ABONDS.DAT
  #---------------------------------
  data <- sprintf("%3s%6.2f", elems, abonds.fin)
  write(data, file = 'abonds.dat', sep = '\n', append = F)
  write(c(1, 1), file = 'abonds.dat', sep = '\n', append = T)
  
  #---------------------------------
  # WRITE DISSOC.DAT
  #---------------------------------
  data2 <- sprintf("%-2s%8d%10.3f%5d%5d%9.2f", elems.dissoc, dissoc.c2, 
                   dissoc.c3, dissoc.c4, dissoc.c5, abonds.dissoc.fin)
  write(data1.dissoc, file = 'dissoc.dat', sep = '\n', append = F)
  write(data2, file = 'dissoc.dat', sep = '\n', append = T)
  write(data3.dissoc, file = 'dissoc.dat', sep = '\n', append = T)
}


##################################
# FUNCTION  CREATE.MODEL
##################################
create.model <- function(Teff_s, logg_s, met_grade_s, logfile){
  #-------------------------------------------------------------
  ################### CREATE FILE MODELES.MOD ##################
  # GRIDS:
  # -0.5 < logg < 1.0, 3300 < Teff < 3900 K, SPHERICAL
  #  1.0 < logg < 3.0, 2600 < Teff < 5250 K, SPHERICAL
  #  3.0 < logg < 5.0, 2900 < Teff < 8000 K, PLANE-PARALLEL
  #-------------------------------------------------------------
  dir = './atm_models/'
  out.of.grid = 0
  
  if(logg_s < 1.0){
    system(paste('cp ', dir, 'innewmarcs2_grid_giants_g_lt_1 innewmarcs2', sep = ''))
    system(paste('cp ', dir, 'marcs2009z-0.50_a+0.20_giants_g_lt_1.mod marcs2009z-0.50.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z-0.25_a+0.10_giants_g_lt_1.mod marcs2009z-0.25.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.00_a+0.00_giants_g_lt_1.mod marcs2009z0.00.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.25_a+0.00_giants_g_lt_1.mod marcs2009z0.25.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.50_a+0.00_giants_g_lt_1.mod marcs2009z0.50.mod', sep = ''))
    if(Teff_s < 3300 | Teff_s > 3900 | met_grade_s < (-0.5)  | met_grade_s > 0.5){
      out.of.grid = 1
    }
  }
  
  if(1.0 <= logg_s & logg_s < 3.0){
    system(paste('cp ', dir, 'innewmarcs2_grid_giants_g_gt_1 innewmarcs2', sep = ''))
    system(paste('cp ', dir, 'marcs2009z-0.50_a+0.20_giants_g_gt_1.mod marcs2009z-0.50.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z-0.25_a+0.10_giants_g_gt_1.mod marcs2009z-0.25.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.00_a+0.00_giants_g_gt_1.mod marcs2009z0.00.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.25_a+0.00_giants_g_gt_1.mod marcs2009z0.25.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.50_a+0.00_giants_g_gt_1.mod marcs2009z0.50.mod', sep = ''))
    if(Teff_s < 2600 | Teff_s > 5250 | met_grade_s < (-0.5) | met_grade_s > 0.5){
      out.of.grid = 1
    }
  }
  
  if(3.0 <= logg_s){
    system(paste('cp ', dir, 'innewmarcs2_grid_dwarfs2 innewmarcs2', sep = ''))
    system(paste('cp ', dir, 'marcs2009z-0.50_a+0.20_dwarfs2.mod marcs2009z-0.50.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z-0.25_a+0.10_dwarfs2.mod marcs2009z-0.25.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.00_a+0.00_dwarfs2.mod marcs2009z0.00.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.25_a+0.00_dwarfs2.mod marcs2009z0.25.mod', sep = ''))
    system(paste('cp ', dir, 'marcs2009z+0.50_a+0.00_dwarfs2.mod marcs2009z0.50.mod', sep = ''))
    if(Teff_s < 2900 | Teff_s > 8000 | met_grade_s < (-0.5)  | met_grade_s > 0.5){
      out.of.grid = 1
    }
  }
  
  if(logg_s < (-0.5) | logg_s > 5.0){
    out.of.grid = 1
  }
  
  if(out.of.grid == 1){
    data <- paste('    ***** WARNING: PARAMETERS OUT OF GRID *****')
    write(data, file = logfile, append = T)
  }
  
  l1 = 'name="modeles"'
  l2 = 'rm -f ${name}.mod'
  l3 = 'rm -f ${name}mod.log'
  l4 = 'touch ${name}.mod'
  l5 = 'touch ${name}mod.log'
  l6 = paste('metl="', met_grade_s,  '"')
  l7 = paste('Teff="', Teff_s,  '"')
  l8 = paste('logg="', logg_s,  '"')
  
  
  l9 = "echo '----------- STARTING THE MEGA-CYCLE -----------'
  lin=0
  lin=`expr $lin + 1`
  ./innewmarcs2 << DONE
  '${name}.mod' '${name}.dat'
  T
  ${Teff} ${logg} ${metl} newmarcs
  ${lin}
  F
  DONE"
 
  data = c(l1, l2, l3, l4, l5, l6, l7, l8, l9)
  
  write(data, file = 'summed.sh', sep = '\n')
  
  system("chmod 777 summed.sh")
  system("./summed.sh")

}

##################################
# FUNCTION SET.STSPEC.FILENAME
##################################
set.stspec.filename <- function(feh, afe, lmin, lmax, Teff, logg, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe){
  temp = c(feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
  temp_s = vector()
  for(j in 1:length(temp)){
    if(temp[j] < 0){
      temp_s[j] = sprintf('%1s%4.2f', '-', abs(temp[j]))
    }else{
      temp_s[j] = sprintf('%1s%4.2f', '+', temp[j])
    }
  }
  
  Teff_s = sprintf('%4.0f', Teff)
  if(logg < 0){
    logg_s = sprintf('%1s%4.2f', '-', abs(logg))
  }else{
    logg_s = sprintf('%1s%4.2f', '+', logg)
  }
  
  lmin_s = sprintf('%4.0f', lmin)
  if(lmax > 1e4){lmax_s = sprintf('%5.0f', lmax)}else{lmax_s = sprintf('%4.0f', lmax)}
  
  file = paste('./Stellar_Spectra/flux_Fe', temp_s[1],  '_a', temp_s[2], '_C', temp_s[3], '_N', temp_s[4],
                      '_O', temp_s[5],'_Mg', temp_s[6], '_Si', temp_s[7],'_Ca', temp_s[8],'_Ti', temp_s[9], 
                     '_Na', temp_s[10], '_T', Teff_s, '_g', logg_s, '_', lmin_s, '-', lmax_s, sep = '')
  return(file)
}

