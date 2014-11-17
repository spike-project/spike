from util.read_msh5 import read_msh5

filemsh5 = '/Users/chiron/FTICR_data/170314_SubP_IRMPD_2D_VD_2D_000001_mr_200314.msh5'
data = read_msh5(filemsh5).resmin
              
def show(kind = 'm/z', scale = 5):
    data.units = kind
    data.display(show = True, scale = scale)

if __name__ == '__main__':
    show()