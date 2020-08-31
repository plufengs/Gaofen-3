# -*- coding: utf-8 -*-
import datetime
from osgeo import gdal,osr
import numpy as np
"""
通用函数
"""
def time_dif(Time1,Time2):
    '''
    calculate the difference between time2 and time1，time:'XXXX-XX-XX XX:XX:XX.XXXXXX'
    :param Time1: 起始时间，
    :param Time2: 结束时间
    :return: 时间差，精确到微秒
    '''
    time1 = datetime.datetime.strptime(Time1,'%Y-%m-%d %H:%M:%S.%f')
    time2 = datetime.datetime.strptime(Time2,'%Y-%m-%d %H:%M:%S.%f')
    return (time2-time1).total_seconds()
'''
对于遥感图像而言，图像行列号（i,j）→投影坐标（不同投影具有不同的坐标表示）→地理坐标（lon,lat）
图像行列号与投影坐标：通过仿射矩阵进行联系
投影坐标与地理坐标：通过投影方式进行联系（比如WGS84地理坐标系进行UTM投影，CGCS2000地理坐标系进行高斯投影）
'''
def geo2lonlat(dataset, x, y):    
    '''
    将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param x: 投影坐标x
    :param y: 投影坐标y
    :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
    '''
    prosrs = osr.SpatialReference()#给定数据的投影参考系
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()#给定数据的地理参考系
    ct = osr.CoordinateTransformation(prosrs, geosrs)
    coords = ct.TransformPoint(x, y)
    return coords[:2]
def lonlat2geo(dataset, lon, lat):
    '''
    将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param lon: 地理坐标lon经度
    :param lat: 地理坐标lat纬度
    :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    prosrs = osr.SpatialReference()#给定数据的投影参考系
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()#给定数据的地理参考系
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lon, lat)
    return coords[:2]
def imagexy2geo(dataset, row, col):
    '''
    根据GDAL的六参数模型将影像图上坐标（行列号）转为投影坐标或地理坐标（根据具体数据的坐标系统转换）
    :param dataset: GDAL地理数据
    :param row: 像素的行号
    :param col: 像素的列号
    :return: 行列号(row, col)对应的投影坐标或地理坐标(x, y)
    '''
    trans = dataset.GetGeoTransform()
    px = trans[0] + col * trans[1] + row * trans[2]
    py = trans[3] + col * trans[4] + row * trans[5]
    return px, py
def geo2imagexy(dataset, x, y):
    '''
    根据GDAL的六参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
    :param dataset: GDAL地理数据
    :param x: 投影或地理坐标x
    :param y: 投影或地理坐标y
    :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col)
    '''
    trans = dataset.GetGeoTransform()
    numerator = trans[1]*trans[5] - trans[2]*trans[4]
    col = int((trans[5]*(x - trans[0]) - trans[2]*(y - trans[3])) / numerator + 0.5)
    row = int((trans[1]*(y - trans[3]) - trans[4]*(x - trans[0])) / numerator + 0.5)
    return (row, col)  
def imagexy2lonlat(dataset,row, col):
    '''
    影像行列转经纬度：
    ：通过影像行列转平面坐标
    ：平面坐标转经纬度
    '''
    coords = imagexy2geo(dataset, row, col)
    coords2 = geo2lonlat(dataset,coords[0], coords[1])
    return (coords2[0], coords2[1])
def lonlat2imagexy(dataset,lon, lat):
    '''
    经纬度转影像行列：
    ：通过经纬度转平面坐标
    ：平面坐标转影像行列
    '''
    coords = lonlat2geo(dataset, lon, lat) 
    coords2 = geo2imagexy(dataset,coords[0], coords[1])
    return (int(round(abs(coords2[0]))), int(round(abs(coords2[1]))))
def clipBylonlat(dataset,lon_min,lon_max,lat_min,lat_max):
    '''
    通过经纬度对影像进行裁剪
    :param dataset: GDAL地理数据
    :param lon_max:最大经度
    :param lon_min: 最小经度
    :param lat_max:最大纬度
    :param lat_min: 最小纬度
    :return: GDAL地理数据
    '''
    row_ul,col_ul = lonlat2imagexy(dataset,lon_min,lat_max)#左上角经纬度对应行列号,即裁剪图像左上角
    row_dr,col_dr = lonlat2imagexy(dataset,lon_max,lat_min)#右下角纬度对应行列号，即裁剪图像右下角
    if col_ul > dataset.RasterXSize:#判断裁剪范围是否超出源图像范围
        col_ul = 0
    if row_ul > dataset.RasterYSize:
        row_ul = 0
    if col_dr > dataset.RasterXSize:
        col_dr = dataset.RasterXSize
    if row_dr > dataset.RasterYSize:
        row_dr = dataset.RasterYSize
    width = col_dr-col_ul#裁剪图像宽度;由于计数是从0开始的，所以不需要+1
    height = row_dr-row_ul#裁剪图像高度
    newDataset = gdal.GetDriverByName(dataset.GetDriver().GetDescription()).Create('file',width,height,dataset.RasterCount,dataset.GetRasterBand(1).DataType)
    GeoT = dataset.GetGeoTransform()
    px = GeoT[0] + col_ul * GeoT[1] + row_ul * GeoT[2]
    py = GeoT[3] + col_ul * GeoT[4] + row_ul * GeoT[5]
    newGeoT = (px,GeoT[1],GeoT[2],py,GeoT[4],GeoT[5])
    newDataset.SetGeoTransform(newGeoT)
    newDataset.SetProjection(dataset.GetProjection())

    if dataset.RasterCount == 1:
        data = dataset.GetRasterBand(1).ReadAsArray(col_ul,row_ul,width,height)
        newDataset.GetRasterBand(1).WriteArray(data)  #写入数组数据
    else:
        for i in range(0,dataset.RasterCount):
            data = dataset.GetRasterBand(i+1).ReadAsArray(col_ul,row_ul,width,height)
            newDataset.GetRasterBand(i+1).WriteArray(data[i])
    return newDataset
def writeTiff(data,filename,description,geotrans,projection):  
    '''
    :param data:遥感影像矩阵
    :param description:遥感影像的描述
    :param geotrans:仿射变换矩阵
    :param projection:投影矩阵
    '''

    #判断写入数据对应的gdal数据类型
    if 'int8' in data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in data.dtype.name:
        datatype = gdal.GDT_UInt16
    elif 'int32' in data.dtype.name:
        datatype = gdal.GDT_UInt32
    elif 'float32' in data.dtype.name:
        datatype = gdal.GDT_Float32
    else:
        datatype = gdal.GDT_Float64
    
    if len(data.shape) == 3:
        bandnums, height, width = data.shape
    else:
        bandnums, (height, width) = 1,data.shape
    
    driver = gdal.GetDriverByName(description)
    dataset = driver.Create(filename, width, height, bandnums, datatype)
    dataset.SetGeoTransform(geotrans) #写入仿射变换参数
    dataset.SetProjection(projection) #写入投影
    
    if bandnums == 1:
        dataset.GetRasterBand(1).WriteArray(data)  #写入数组数据
    else:
        for i in range(bandnums):
            dataset.GetRasterBand(i+1).WriteArray(data[i])
    del dataset
def linear5(imArr):
    '''
    影像5%线性拉伸
    :param imArr:输入影像矩阵
    :return outArr:输出影像矩阵
    '''
    #5%线性拉伸
    tmp = imArr
    tmp = tmp.flatten()
    tmp = np.delete(tmp,np.where(tmp==0.0))
    p1 = np.percentile(tmp,5)
    p2 = np.percentile(tmp,95)
    maxout=255
    minout=0    
    outArr=minout + ( (imArr-p1) / (p2-p1) ) * (maxout - minout)
    outArr[outArr < minout]=minout
    outArr[outArr > maxout]=maxout
    outArr[np.where(imArr==0.0)]=maxout
    return outArr