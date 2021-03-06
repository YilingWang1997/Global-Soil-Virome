import gdal
import numpy as np
import sys
import pandas as pd

data=np.loadtxt("cal_index.csv",delimiter=",") #Input the matrixs generated by R

file1=gdal.Open('$factor.a.tif')
band1=file1.GetRasterBand(1)
factor_a=band1.ReadAsArray()

file1=gdal.Open('$factor.b.tif')
band1=file1.GetRasterBand(1)
factor_b=band1.ReadAsArray()    ##Input layer

factor_a[1,1] #The number in this position is likely to be "NA".
#If the output is "NA", skip the two steps below.
factor_a[factor_a==factor_a[1,1]]=0 
factor_a[factor_a!=0]=1  

factor_b[1,1] 
factor_b[factor_b==factor_b[1,1]]=0
factor_b[factor_b!=0]=1

file1=gdal.Open('New_TotalNumber_Bootstrap_CoefVar.tif')
band1=file1.GetRasterBand(1)
maskocean = band1.ReadAsArray() 
maskocean[np.isnan(maskocean)] = maskocean[1,1]
maskocean[maskocean==maskocean[1,1]]=0
maskocean[maskocean!=0]=1

mask=factor_a*factor_b*maskocean

data_mask=data*mask

#Map generation
path1='$input.tif'#choose one tif with customized coordinate system
output_path=#Reference:"index.tif"
driver=gdal.GetDriverByName('GTiff')
file1=gdal.Open(path1)
geotrans=file1.GetGeoTransform()
proj=file1.GetProjection()
band1=file1.GetRasterBand(1)
file2=driver.Create(output_path,band1.XSize,band1.YSize,1,band1.DataType)
file2.SetGeoTransform(geotrans)
file2.SetProjection(proj)
band2=file2.GetRasterBand(1)
Array_New=data_mask
band2.WriteArray(Array_New)
file2=None