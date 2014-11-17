from util.debug_tools import*  
import File.Apex as Apex
import numpy as np
from NPKData import NPKData, as_cpx

@decclassdebugging                               
class FTICR_treat():
    def __init__(self):
        pass

    def read(self, add):
        '''
        Reads Apex or msh5 files. If Apex is used makes locally a corresponding msh5.
        '''
        print "path exists? ",os.path.exists(add)
        nameroot = os.path.basename(add).split(".")[0]
        splitadd = os.path.splitext(add)
        ext = splitadd[-1]
 
        if ext == '.d' and not os.path.exists(splitadd[0]+".msh5"):
            data = Apex.Import_1D(add, '')
            print "just after Apex kind of data is ", type(data)
        elif ext == '.msh5' or os.path.exists(splitadd[0]+".msh5"):
            filemsh5 = nameroot+'.msh5'
        else : 
            print "Unknown format"

        return data
    
    def write(self, data, filename):
        print "writing to hdf5 file "
        print 'filename is ', filename
        print "dir(data) ", dir(data)
        #data.save(filename)
        data.save_msh5(filename)

    def xymz(self, datafticr, specfreq, mzmax, resol = 1):
        '''
        transformation for plotting in m/z
        specmzindex  contains the m/z values, x axis
        specmz the corresponding values, y axis
        datafticr (DataFTICR format) is used for recuperating itomz
        '''
        print "in xymz "
        specmz = specfreq
        print "specfreq.size ",specfreq.size
        specmzindex = np.zeros(specfreq.size-1) # initialization
        if resol != 1 : 
            print "superresolution"
        for i in range(1, specmzindex.size): # Super resolution size in case of super resol.. 
            mzval = datafticr.axes(1).itomz(i)
            if resol != 1 : 
                mzcorrect = datafticr.axes(1).itomz(i/float(resol))
                specmzindex[i-1] =  mzcorrect
            else :
                specmzindex[i-1] =  mzval
        print "specmzindex[specmzindex < mzmax]",specmzindex[specmzindex < mzmax]
        indexmz = np.where(specmzindex < mzmax)[0]
        yspecmz  = specmz[1:][indexmz]
        xmz = specmzindex[indexmz]
        yspecmz[np.where(xmz < 300)[0]] = 0 
        return xmz, yspecmz

@decclassdebugging                               
class Orbitrap_treat():
    def __init__(self):
        self.refnorm = 1
        self.refresol = 1
    
    def find_range(self,f):
        '''
        Calculate factor self.range_fact in function of the masss range.
        '''
        a = os.path.basename(f)
        try:
            result = int(re.search('range(\d*)-.*', a).group(1))
            print "result is ", result
            if result <= 100 :
                self.range_fact = 1
            else:
                self.range_fact = 0.5 
        except Exception:
            self.range_fact = 0.5    

    def read(self, f):
        """
        reads Orbitrap .dat FID files
        returns a numpy array with the FID
        """
        print "f ", f
        self.find_range(f)
        with open(f,'rb') as F:
            self.entete = ''
            pos = 0
            for l in F:
                self.entete += l
                pos += len(l)
                if l.startswith("Data Points"):
                    print re.findall(r'\d+', l)[0]
                if l.startswith("Data:"):
                    break
            F.seek(pos)
            data_interm = np.array(np.fromfile(F, dtype = 'f4').tolist())  # ndarray to list to np.array
            data = NPKData(buffer = data_interm )
        return data
    
    def write(self, data, filename):
        f = open(filename,'wb')
        print "in write orbitrap type(data) is ", type(data)
        new_data = self.entete+data.astype('f4').tostring()
        f.write(new_data)
        f.close()
         
    def mzshow(self, spec, col = 'k', linestyle ='-'):
        '''
        print spectrum in m/z
        '''
        print "spec.size",spec.size
        point_ref = spec.size*377500./(2*999982)
        print "point_ref ",point_ref
        mz_ref = 715.3122
        trunc = 0.1 * spec.size
        xaxis = mz_ref/(((1 + np.arange(1.0 * spec.size))/point_ref)**2)
        plt.plot(xaxis[trunc:], spec[trunc:], color = col, linestyle = linestyle)
        
    def zerfapod(self, data, zerofill = False):
        '''
        Zerofilling and apodisation
        '''
        print "self.refresol ",self.refresol
        #apod = apod_sq_sin
        #apod = apod_sin
        if zerofill :
            spec = data.apod_sq_sin(maxi = 0.5).chsize(self.refresol).rfft().modulus().buffer
        else :
            spec = data.apod_sq_sin(maxi = 0.5).rfft().modulus().buffer
        print "spec.size",spec.size
        return spec
    
    def xymz(self, y):
        '''
        Prepare data for orbitrap m/z unit.
        '''
        mz_ref = 715.3122 # reference mass
        trunc_deb = 0.1
        data_size = y.size
        point_ref = data_size * 377500./(2*999982)
        trunc = trunc_deb * data_size
        xaxis = mz_ref/(((1 + np.arange(1.0*data_size)*self.range_fact)/point_ref)**2)
        xmz = xaxis[trunc:]
        yspecmz =  y[trunc:]
        print "xmz, yspecmz ", xmz, yspecmz
        return xmz, yspecmz