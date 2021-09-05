import gdal
import numpy as np
import gdalconst
import numpy.ma as ma
import os,sys
#for i in range(1,2):

path = '$'#The workdir of layers needed to be resample
filepath = os.listdir(path)

for f in filepath:
    if f[-3:] == 'tif':
        print('%sStart Calculating'%f)
        path2 = path + f
        i = 1
        if i<10:
            path_New = '$/New_' + '%s'%f
        else:
            path_New = '$'  + '/New_'+ path[-20:-4] + ".tif"
        file1 = gdal.Open(path2)
        geotrans = file1.GetGeoTransform()
        band1 = file1.GetRasterBand(1)
        print('仿射矩阵如下')
        print(geotrans)

        inputProj = file1.GetProjection()
        referencefileProj = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]' 
        outbounding = [-180,-90,180,90]
        #file2 = gdal.Warp(path_New,path2,xRes = 0.01,yRes = 0.01,creationOptions = ['TILED = YES'],outputBounds=outbounding,srcSRS=inputProj, dstSRS=referencefileProj,outputType = gdalconst.GDT_Float32)
        file2 = gdal.Warp(path_New,path2,xRes = 0.01,yRes = 0.01,creationOptions = ['TILED = YES'],outputBounds=outbounding, resampleAlg=gdalconst.GRA_Bilinear,srcSRS=inputProj, dstSRS=referencefileProj,outputType = gdalconst.GDT_Float32)
        #path_New = None
        file1 = None
        band1 = None

        #file2 = gdal.Open(path_New,1)
        band2 = file2.GetRasterBand(1)
        p = 0#
        #p = float(band2.GetNoDataValue())
        print(p)
        Array2 = band2.ReadAsArray() 
        if (p<0):
            p = np.min(Array2)
        elif (p==0):
            p = 0
        elif (p==65536):
            p = p
        else:
            p = np.max(Array2)
        print(p)
        #Array2 = Array2.astype('float64')
        mask = Array2 == p
        Array3 = ma.masked_array(Array2,mask=mask)
        Mean = np.nanmean(Array3)
        print('Mean：' + '%s'%Mean)
        Std = np.nanstd(Array2)
        print('Std：' + '%s'%Std)
        print('Original max：' + '%s'%np.nanmax(Array2))
        print('Original min：' + '%s'%np.nanmin(Array2))
        Array_New = (Array2 - Mean)/Std
        print('Max：' + '%s'%np.nanmax(Array_New))
        print('Min：' + '%s'%np.nanmin(Array_New))
        Array_New = Array_New.astype('float32')
        band2.WriteArray(Array_New)
        print('Graph size：')
        print(Array_New.shape)
        file2 = None
        print('%sDone'%f)
'''
for t in range(1):
    path = '/data1/ylwang/Gisdata/soilbiomass/c_1km/'
    filepath = os.listdir(path)
    for f in filepath:
        if f == 'hdr.adf':
            print('bio_%sStart Calculating'%t)
            path2 = path + f
            i = 1
            if i<10:
                path_New = '/data1/ylwang/Gisdata/TC_New/New_soilbiomass_1km.tif'
            else:
                path_New = '/data1/ylwang/Gisdata/TC_New'  + '/New_'+ path[-20:-4] + ".tif"
            file1 = gdal.Open(path2)
            geotrans = file1.GetGeoTransform()
            band1 = file1.GetRasterBand(1)

            inputProj = file1.GetProjection()
            referencefileProj = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]' 
            outbounding = [-180,-90,180,90]
            file2 = gdal.Warp(path_New,path2,xRes = 0.01,yRes = 0.01,outputBounds=outbounding, resampleAlg=gdalconst.GRA_Bilinear,creationOptions = ["COMPRESS = LZW"],srcSRS=inputProj, dstSRS=referencefileProj,outputType = gdalconst.GDT_Float32,workingType = gdalconst.GDT_Float64)
            #path_New = None
            file1 = None
            band1 = None

            #file2 = gdal.Open(path_New,1)
            band2 = file2.GetRasterBand(1)
            p = float(band2.GetNoDataValue())
            print(p)
            Array2 = band2.ReadAsArray() 
            if (p<0):
                p = np.min(Array2)
            elif (p==0):
                p = 0
            elif (p==65536):
                p = p
            else:
                p = np.max(Array2)
            print(p)
            Array2 = Array2.astype('float64')
            mask = Array2 == p
            Array2 = ma.masked_array(Array2,mask=mask)
            Mean = np.nanmean(Array2)
            print('Mean：' + '%s'%Mean)
            Std = np.nanstd(Array2)
            print('Std：' + '%s'%Std)
            print('Original max：' + '%s'%np.nanmax(Array2))
            print('Original min：' + '%s'%np.nanmin(Array2))
            Array_New = (Array2 - Mean)/Std
            print('Max：' + '%s'%np.nanmax(Array_New))
            print('Min：' + '%s'%np.nanmin(Array_New))
            Array_New = Array_New.astype('float32')
            band2.WriteArray(Array_New)
            print('Graph size：')
            print(Array_New.shape)
            print('%sDone'%f)
'''

