import astropy.io.fits as pyfits, numpy as np, sys, os, pickle
DIR='/archive_data/desarchive/DEC/finalcut/Y5A1/HITS/3505/'

def listdir_fullpath(d):
  return [os.path.join(d, f) for f in os.listdir(d)]
def distance(ra1,dec1,ra2,dec2):
  return np.sqrt((np.cos(np.pi/360.*(dec1+dec2))*(ra1-ra2))**2+(dec1-dec2)**2)
def opensurvey(inname):
  return pickle.load(open(inname, "rb" ) )
def getcoaddindex(ras,decs,rac,decc,r=1.5/3600.):
  cosdecs = np.cos(np.pi*decs/180.)
  cosdecs2 = cosdecs**2
  ns = decs.size
  nc = decc.size
  coadd_index = np.ones(ns)*-1
  if (nc == 0) or (ns == 0):
    return coadd_index.astype(int)
  ncmin = 0
  for num in range(ns):
    [ra,dec,cosdec,cosdec2] = [ras[num], decs[num], cosdecs[num], cosdecs2[num]]
    decmin = dec-r

    while (decc[ncmin] < decmin) & (ncmin < nc-1):
      ncmin+=1
    num2 = ncmin
    rmin2 = r**2
    while (decc[num2] < dec+r):
      ddec2 = (decc[num2]-dec)**2
      if ddec2 < rmin2:
        r2 = ddec2+(rac[num2]-ra)**2*cosdec2
        if r2 < rmin2:
          rmin2 = r2
          coadd_index[num] = num2
      num2+=1
      if num2 == nc:
        break
  return coadd_index.astype(int)
  

class image:
  def __init__(self,EXPDIR):
    self.expdir = EXPDIR
    cats = np.sort(listdir_fullpath(self.expdir))
    self.cats=cats[np.array(['fullcat' in cat for cat in cats])]
  def stats(self):
    print self.expdir
    self.minras=[]
    self.maxras=[]
    self.mindecs=[]
    self.maxdecs=[]
    self.nsources=[]
    for cat in self.cats:
      table=pyfits.open(cat)[2]  
      RA=table.data['ALPHA_J2000']
      DEC=table.data['DELTA_J2000']
      self.minras.append(np.min(RA[RA > 0])) 
      self.maxras.append(np.max(RA[RA < 360]))
      self.mindecs.append(np.min(DEC[DEC > -90]))
      self.maxdecs.append(np.max(DEC[DEC < 90])) 
      self.nsources.append(table.data.shape[0])
    self.minras = np.array(self.minras)
    self.maxras = np.array(self.minras)
    self.mindecs = np.array(self.mindecs)
    self.maxdecs = np.array(self.mindecs)
    self.nsources = np.array(self.nsources)
    if self.minras.size > 0:
      self.band=pyfits.open(self.cats[0])[0].header['BAND']
      self.minra = np.min(self.minras)
      self.maxra = np.max(self.maxras)
      self.mindec = np.min(self.mindecs)
      self.maxdec = np.max(self.maxdecs)
      self.nsource = np.sum(self.nsources)
      self.ra = np.mean(self.minra,self.maxra)    
      self.dec = np.mean(self.mindec,self.maxdec)    
    else:
      self.band='None'
      self.minra = np.nan
      self.maxra = np.nan
      self.mindec = np.nan
      self.maxdec = np.nan
      self.nsource = 0
      self.ra = np.nan 
      self.dec = np.nan 

class field:
  def __init__(self,EXPDIRS,DIR):
    self.images=[]
    for EXPDIR in EXPDIRS:
      self.images.append(image(EXPDIR))
    self.nimages = len(self.images)
    self.outdir = DIR
  def cat_combine(self,image, outname):
    if os.path.exists(outname):
      print outname+' already exists. Not remaking.'
      return pickle.load(open(outname,'rb'))
    [flags, ra, dec, ra_err, dec_err, mag, mag_err, class_star] = [[],[],[],[],[],[],[],[]]
    for cat in image.cats:
      fits= pyfits.open(cat)[2].data
      flags.append(fits['FLAGS']); ra.append(fits['ALPHA_J2000']); dec.append(fits['DELTA_J2000']); ra_err.append(fits['ERRAWIN_WORLD']); dec_err.append(fits['ERRBWIN_WORLD']); mag.append(fits['MAG_PSF']); mag_err.append(fits['MAGERR_PSF']); class_star.append(fits['CLASS_STAR'])
    flags=np.concatenate(flags); ra=np.concatenate(ra); dec=np.concatenate(dec); ra_err=np.concatenate(ra_err); dec_err=np.concatenate(dec_err); mag=np.concatenate(mag); mag_err=np.concatenate(mag_err); class_star=np.concatenate(class_star)
    good = (class_star > 0.8) & (flags == 0) & (mag_err < 0.1) & (mag_err > 0.01)
    ra=ra[good]; dec = dec[good]; ra_err=ra_err[good]; dec_err=dec_err[good]; mag = mag[good]; mag_err = mag_err[good]
    order = np.argsort(dec)
    ra=ra[order]; dec = dec[order]; ra_err=ra_err[order]; dec_err=dec_err[order]; mag = mag[order]; mag_err = mag_err[order]
    output=[ra,dec,ra_err,dec_err,mag,mag_err]
    print outname+' written.'
    pickle.dump(output,open(outname,'wb'))
    return output
  def cat_combines(self):
    self.stars=[]
    for image in self.images:
      self.stars.append(self.cat_combine(image,self.outdir+'/'+image.expdir.split('/')[-3]+'.pkl'))
  def crossmatch(self):
    outname = self.outdir+'/crossmatch.pkl'
    if os.path.exists(outname):
      print outname+' already exists. Not remaking.'
      return pickle.load(open(outname,'rb'))
    self.nstars = []
    for star in self.stars:
      self.nstars.append(len(star[0]))
    self.nstars = np.array(self.nstars)
    [self.nstars_max, self.argmax] = [np.max(self.nstars), np.argmax(self.nstars)]
    [ra, dec, ra_err, dec_err, mag, mag_err] = [np.zeros([len(self.stars), self.nstars_max]), np.zeros([len(self.stars), self.nstars_max]), np.zeros([len(self.stars), self.nstars_max]), np.zeros([len(self.stars), self.nstars_max]), np.zeros([len(self.stars), self.nstars_max]), np.zeros([len(self.stars), self.nstars_max])]
    for num in range(len(self.stars)):
      if num == self.argmax:
        [ra[num], dec[num], ra_err[num], dec_err[num], mag[num], mag_err[num]] = [self.stars[num][0], self.stars[num][1], self.stars[num][2], self.stars[num][3], self.stars[num][4], self.stars[num][5]]
        continue
      index = getcoaddindex(self.stars[num][0], self.stars[num][1], self.stars[self.argmax][0], self.stars[self.argmax][1]) 
      matched = (index > -1)   
      ra[num][index[matched]] = self.stars[num][0][matched] 
      dec[num][index[matched]] = self.stars[num][1][matched] 
      ra_err[num][index[matched]] = self.stars[num][2][matched] 
      dec_err[num][index[matched]] = self.stars[num][3][matched] 
      mag[num][index[matched]] = self.stars[num][4][matched] 
      mag_err[num][index[matched]] = self.stars[num][5][matched] 
    output = [ra,dec,ra_err,dec_err,mag, mag_err]
    pickle.dump(output,open(outname,'wb'))
    return output
  def makeplots(self):
    cm = self.crossmatch()
    diff = cm[4] - cm[4][11]
    good = (cm[4] > 0) & (diff != 0)
    x = cm[5][good]
    y = diff[good]

class survey:
  def __init__(self,DIR):
    EXPDIRS=[]
    for DIR1 in listdir_fullpath(DIR):
      for DIR2 in listdir_fullpath(DIR1):
        if os.path.exists(DIR2+'/p01/cat'):
          EXPDIRS.append(DIR2+'/p01/cat')
    self.expdirs = EXPDIRS
    self.images=[]
    for EXPDIR in EXPDIRS:
      self.images.append(image(EXPDIR))
    self.calc_stats = False
    self.nimages = len(self.images)
  def stats(self):
    num =0 
    self.bands=[]
    for image in self.images:
      num += 1
      print str(num)+' of '+str(self.nimages)
      image.stats()
      self.bands.append(image.band)
    self.bands = np.array(self.bands)
    self.calc_stats = True
  def save(self, outname):
    pickle.dump(self,open(outname,'wb'))
  def fields(self, minrad = 0.5):
    self.fieldras=[]
    self.fielddecs=[]
    self.nfields=0
    self.fieldids = -1*np.ones(self.nimages)
    for num in range(self.nimages):
      image = self.images[num]
      if image.nsource == 0:
        continue
      if self.nfields==0:
        self.fieldras  = np.append(self.fieldras,  image.ra)
        self.fielddecs = np.append(self.fielddecs, image.dec)
        self.fieldids[num] = self.nfields
        self.nfields += 1
        continue
      distances = distance(self.fieldras, self.fielddecs, image.ra, image.dec) 
      if np.min(distances) < minrad:
        self.fieldids[num] = np.argmin(distances)
      else:
        self.fieldras  = np.append(self.fieldras,  image.ra)
        self.fielddecs = np.append(self.fielddecs, image.dec)
        self.fieldids[num] = self.nfields
        self.nfields += 1
  def makdirs(self, ROOTDIR='./'):
    for field in np.unique(self.fieldids):
      DIR1 = ROOTDIR+"/F%d" %field
      os.mkdir(DIR1)
      for band in np.unique(self.bands[self.fieldids == field]):
        if band != 'None':
          DIR2 = DIR1+'/'+band
          print DIR2, np.sum((self.bands == band) & (self.fieldids == field))
          os.mkdir(DIR2)
