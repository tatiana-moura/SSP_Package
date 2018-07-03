#########################
# TEST SCRIPT
#########################
source('stpars.R')
source('pfant12.R')
source('SSP_model.R')

n_ms = 2
n_rg = 2
feh = -1.0
afe = 0.40
age = 10
lmin = 4601.60
lmax = 5600.60
imf = 'Kroupa'
dl = 0.1

stpars(n_ms, n_rg, feh, afe, age, fig = T)
pfant12(feh, afe, lmin, lmax, age, n_ms = n_ms, n_rg = n_rg, dl = dl)
ssp.model(feh, afe, age, imf, fwhm = 2.5, n_ms = n_ms, n_rg = n_rg, 
          lmin = lmin, lmax = lmax, slope = slope, dl = dl)


file_ssp <- set.ssp.filename(feh, afe, age = age, imf = imf, slope = slope)
t = read.table(file_ssp)
xx1 = median(t$V2[t$V1 > lmin & t$V1 < lmax])
plot(t$V1, t$V2, col = 'black', type = 'l', lty = 1, ylab = 'Flux', 
     xlab = expression(paste(lambda, ' (', ring(A), ')')))

t = read.table('Iku1.30Zp0.22T08.9125')
xx2 = median(t$V2[t$V1 > lmin & t$V1 < lmax])
lines(t$V1, t$V2/(xx2/xx1), type = 'l', col ='red', lwd = 0.8)

legend(8400, 5e-5, c('Synthetic SSP', 'MIUSCAT'), bty = 'n', lty = c(1, 1),
       lwd = c(1, 0.8), col = c('black', 'red'))




