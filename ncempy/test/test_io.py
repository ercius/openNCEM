import ncempy.io as nio
from pathlib import Path

dPath = Path(r'C:\Users\linol\Data')

print('ser')
fPath = Path('01_Si110_5images_1.ser')
_ = nio.ser.serReader(dPath / fPath)
print(_.keys())

print('dm')
fPath = Path('5_Te_15x83_ss=3nm_CL=245_alpha=p48_p06sec_no beamstop_bin4_300kV.dm4')
_ = nio.dm.dmReader(dPath / fPath)
print(_.keys())

print('emd')
fPath = Path('TimeSeries_20.emd')
_ = nio.emd.emdReader(dPath / fPath)
print(_.keys())

print('emdVelox')
fPath = Path('1435 1.2 Mx STEM HAADF-DF4-DF2-BF.emd')
_ = nio.emdVelox.emdVeloxReader(dPath / fPath)
print(_.keys())

print('mrc')
fPath = Path('AgNWweld_tomo2_1wire_115kx_160mmCL.mrc')
_ = nio.mrc.mrcReader(dPath / fPath)
print(_.keys())
