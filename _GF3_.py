# -*- coding: utf-8 -*-
"""
Class GF3
从GF3头文件中提取所有有关信息
"""

import xml.dom.minidom
from osgeo import gdal
import numpy as np
import _op_
import glob,time,os
'''
通用类
'''
class Sensor:
    def __init__(self):
        self.lamda = 0.0
        self.RadarCenterFrequency= 0.0
class Wave:
    def __init__(self):
        self.centerLookAngle = 0.0
        self.prf = 0.0
        self.sampleRate = 0.0
        self.sampleDelay = 0.0
        self.groundVelocity = 0.0
class Platform:
    def __init__(self):
        self.CenterTime = ''
        self.Rs = 0.0
        self.satVelocity = 0.0
        self.RollAngle = 0.0
        self.PitchAngle = 0.0
        self.YawAngle = 0.0
        self.Xs = 0.0
        self.Ys = 0.0
        self.Zs = 0.0
        self.Vxs = 0.0
        self.Vys = 0.0
        self.Vzs = 0.0    
class GPSParam:
    def __init__(self):
        self.TimeStamp = ''
        self.xPosition = 0.0
        self.yPosition = 0.0
        self.zPosition = 0.0
        self.xVelocity = 0.0
        self.yVelocity = 0.0
        self.zVelocity = 0.0
        self.ZeroTime = 0.0
class ATTIParam:
    def __init__(self):
        self.TimeStamp = ''
        self.yawAngle = 0.0
        self.rollAngle = 0.0
        self.pitchAngle = 0.0
class Productinfo:
    def __init__(self):
        self.NominalResolution = 0.0
        self.WidthInMeters = 0.0 
        self.productLevel = 0
        self.productType = ''
        self.productFormat = ''
        self.productGentime = ''
        self.productPolar = ''
class Imageinfo:
    def __init__(self):
        self.imagingTime_start = ''
        self.imagingTime_end = ''
        self.TimeTotal = 0.0#总的成像时间
        self.TimeInterval = 0.0#每行时间间隔
        self.FirstTime = 0.0#以图像开始成像时间为0时刻
        self.LastTime = 0.0#图像结束时刻，0时刻为图像开始成像时间
        self.nearRange = 0.0
        self.refRange = 0.0
        self.eqvFs = 0.0
        self.eqvPRF = 0.0
        self.center = []#[latitude,longitude]
        self.topLeft = []#~
        self.topRight = []#~
        self.bottomLeft = []#~
        self.bottomRight = []#~
        self.width = 0
        self.height = 0
        self.widthspace = 0.0
        self.heightspace = 0.0
        self.sceneShift = 0.0
        self.imagebit = ''
        self.QualifyValue = []#[HH,HV,VH,VV]
        self.echoSaturation = []#~
class Processinfo:
    def __init__(self):
        self.EphemerisData = ''
        self.AttitudeData = ''
        self.algorithm = ''
        self.Calibration = []#[HH,HV,VH,VV]
        self.AzFdc0 = 0.0
        self.AzFdc1 = 0.0
        self.MultilookRange = 0
        self.MultilookAzimuth = 0
        self.DoIQComp = 0
        self.DoChaComp = 0
        self.RangeWeightType = 0
        self.RangeWeightPara = []#[float,float]
        self.AzimuthWeightType = 0
        self.AzimuthWeightPara = []#[float,float]
        self.EarthModel = ''
        self.ProjectModel = ''
        self.incidenceAngleNearRange = 0.0
        self.incidenceAngleFarRange = 0.0
        self.RangeLookBandWidth = 0.0
        self.AzimuthLookBandWidth = 0.0
        self.TotalProcessedAzimuthBandWidth = 0.0
        self.DopplerParametersReferenceTime = 0.0
        self.DopplerCentroidCoefficients = []#[d0,d1,d2,d3,d4]
        self.DopplerRateValuesCoefficients = []#[r0,r1,r2,r3,r4]
class Incidence:
    def __init__(self):
        self.nums = 0
        self.stepSize = 0
        self.Value = []

class GF3:
    def __init__(self):
        self.Sensor = Sensor()
        self.Wave = Wave()
        self.lookDirection = ''
        self.antennaMode = ''
        self.polarMode = []
        self.polarNums = 0
        self.Platform = Platform()
        self.GPSParam = []#[GPStime,px,py,pz,vx,vy,vz,以图像成像时间为0时刻的GPS时间]
        self.ATTIParam = []
        self.Productinfo = Productinfo()
        self.Imageinfo = Imageinfo()
        self.Processinfo = Processinfo()
        self.Incidence = Incidence()
        self.SatEveryRow = []#每一行对应的卫星状态
        self.Dataset = []
    def read_meta(self,meta_XML):  
        meta_Type = meta_XML.split('_')
        if 'L1A' in meta_Type:
            self.read_L1A(meta_XML)
        else:
            print('Error: Not L1A.XML')
        self.get_satAllRow()
        return 0
    def read_L1A(self,meta_XML):
        meta_dom = xml.dom.minidom.parse(meta_XML)
        root = meta_dom.documentElement

        self.Sensor.lamda = float(root.getElementsByTagName('lamda')[0].childNodes[0].data)
        self.Sensor.RadarCenterFrequency = float(root.getElementsByTagName('RadarCenterFrequency')[0].childNodes[0].data)
        
        self.Wave.centerLookAngle = float(root.getElementsByTagName('centerLookAngle')[0].childNodes[0].data)
        self.Wave.prf = float(root.getElementsByTagName('prf')[0].childNodes[0].data)
        self.Wave.sampleRate = float(root.getElementsByTagName('sampleRate')[0].childNodes[0].data)
        self.Wave.sampleDelay = float(root.getElementsByTagName('sampleDelay')[0].childNodes[0].data)
        self.Wave.groundVelocity = float(root.getElementsByTagName('groundVelocity')[0].childNodes[0].data)
        
        self.lookDirection = root.getElementsByTagName('lookDirection')[0].childNodes[0].data
        self.antennaMode= root.getElementsByTagName('antennaMode')[0].childNodes[0].data
        _polarMode = root.getElementsByTagName('polarMode')[0].childNodes[0].data
        if _polarMode == 'AHV':
            self.polarNums = 4
            self.polarMode = ['HH','HV','VH','VV']
        else:
            self.polarNums = int(_polarMode.__len__()/2)
            for i in range(0,self.polarNums):
                self.polarMode.append(_polarMode[i*2:i*2+2])
            
        self.Platform.CenterTime = root.getElementsByTagName('CenterTime')[0].childNodes[0].data
        self.Platform.Rs = float(root.getElementsByTagName('Rs')[0].childNodes[0].data)
        self.Platform.satVelocity = float(root.getElementsByTagName('satVelocity')[0].childNodes[0].data)
        self.Platform.RollAngle = float(root.getElementsByTagName('RollAngle')[0].childNodes[0].data)
        self.Platform.PitchAngle = float(root.getElementsByTagName('PitchAngle')[0].childNodes[0].data)
        self.Platform.YawAngle = float(root.getElementsByTagName('YawAngle')[0].childNodes[0].data)
        self.Platform.Xs = float(root.getElementsByTagName('Xs')[0].childNodes[0].data)
        self.Platform.Ys = float(root.getElementsByTagName('Ys')[0].childNodes[0].data)
        self.Platform.Zs = float(root.getElementsByTagName('Zs')[0].childNodes[0].data)
        self.Platform.Vxs = float(root.getElementsByTagName('Vxs')[0].childNodes[0].data)
        self.Platform.Vys = float(root.getElementsByTagName('Vys')[0].childNodes[0].data)
        self.Platform.Vzs = float(root.getElementsByTagName('Vzs')[0].childNodes[0].data)
        
        self.Imageinfo.imagingTime_start = root.getElementsByTagName('imagingTime')[0].getElementsByTagName('start')[0].childNodes[0].data
        self.Imageinfo.imagingTime_end = root.getElementsByTagName('imagingTime')[0].getElementsByTagName('end')[0].childNodes[0].data
        self.Imageinfo.nearRange = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('nearRange')[0].childNodes[0].data)
        self.Imageinfo.refRange = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('refRange')[0].childNodes[0].data)
        self.Imageinfo.eqvFs = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('eqvFs')[0].childNodes[0].data)
        self.Imageinfo.eqvPRF = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('eqvPRF')[0].childNodes[0].data)
        self.Imageinfo.center = [float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('center')[0].getElementsByTagName('longitude')[0].childNodes[0].data),\
        float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('center')[0].getElementsByTagName('latitude')[0].childNodes[0].data)]
        self.Imageinfo.topLeft = [float(root.getElementsByTagName('topLeft')[0].getElementsByTagName('longitude')[0].childNodes[0].data),\
        float(root.getElementsByTagName('topLeft')[0].getElementsByTagName('latitude')[0].childNodes[0].data)]
        self.Imageinfo.topRight = [float(root.getElementsByTagName('topRight')[0].getElementsByTagName('longitude')[0].childNodes[0].data),\
        float(root.getElementsByTagName('topRight')[0].getElementsByTagName('latitude')[0].childNodes[0].data)]
        self.Imageinfo.bottomLeft = [float(root.getElementsByTagName('bottomLeft')[0].getElementsByTagName('longitude')[0].childNodes[0].data),\
        float(root.getElementsByTagName('bottomLeft')[0].getElementsByTagName('latitude')[0].childNodes[0].data)]
        self.Imageinfo.bottomRight = [float(root.getElementsByTagName('bottomRight')[0].getElementsByTagName('longitude')[0].childNodes[0].data),\
        float(root.getElementsByTagName('bottomRight')[0].getElementsByTagName('latitude')[0].childNodes[0].data)]
        self.Imageinfo.width = int(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('width')[0].childNodes[0].data)
        self.Imageinfo.height = int(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('height')[0].childNodes[0].data)
        self.Imageinfo.widthspace = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('widthspace')[0].childNodes[0].data)
        self.Imageinfo.heightspace = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('heightspace')[0].childNodes[0].data)
        self.Imageinfo.sceneShift = float(root.getElementsByTagName('imageinfo')[0].getElementsByTagName('sceneShift')[0].childNodes[0].data)
        self.Imageinfo.imagebit = root.getElementsByTagName('imageinfo')[0].getElementsByTagName('imagebit')[0].childNodes[0].data
        _QualifyValue = root.getElementsByTagName('imageinfo')[0].getElementsByTagName('QualifyValue')[0]
        self.Imageinfo.QualifyValue = [_QualifyValue.getElementsByTagName('HH')[0].childNodes[0].data,\
                                       _QualifyValue.getElementsByTagName('HV')[0].childNodes[0].data,\
                                       _QualifyValue.getElementsByTagName('VH')[0].childNodes[0].data,\
                                       _QualifyValue.getElementsByTagName('VV')[0].childNodes[0].data]
        for i in range(0,self.Imageinfo.QualifyValue.__len__()):
            if self.Imageinfo.QualifyValue[i] == 'NULL':
                self.Imageinfo.QualifyValue[i] = 0
            else:
                self.Imageinfo.QualifyValue[i] = float(self.Imageinfo.QualifyValue[i])
        _echoSaturation = root.getElementsByTagName('imageinfo')[0].getElementsByTagName('echoSaturation')[0]
        self.Imageinfo.echoSaturation = [_echoSaturation.getElementsByTagName('HH')[0].childNodes[0].data,\
                                         _echoSaturation.getElementsByTagName('HV')[0].childNodes[0].data,\
                                         _echoSaturation.getElementsByTagName('VH')[0].childNodes[0].data,\
                                         _echoSaturation.getElementsByTagName('VV')[0].childNodes[0].data]
        for i in range(0,self.Imageinfo.echoSaturation.__len__()):
            if self.Imageinfo.echoSaturation[i] == 'NULL':
                self.Imageinfo.echoSaturation[i] = 0
            else:
                self.Imageinfo.echoSaturation[i] = float(self.Imageinfo.echoSaturation[i])
        self.Imageinfo.TimeTotal = _op_.time_dif(self.Imageinfo.imagingTime_start,self.Imageinfo.imagingTime_end)
        self.Imageinfo.TimeInterval = self.Imageinfo.TimeTotal/(self.Imageinfo.height-1)
        self.Imageinfo.LastTime = self.Imageinfo.TimeTotal
        
        _GPSParam = root.getElementsByTagName('GPSParam')
        for i in range(0,_GPSParam.__len__()):
            tmpGPSP = GPSParam()
            tmpGPSP.TimeStamp = _GPSParam[i].getElementsByTagName('TimeStamp')[0].childNodes[0].data
            tmpGPSP.xPosition = float(_GPSParam[i].getElementsByTagName('xPosition')[0].childNodes[0].data)
            tmpGPSP.yPosition = float(_GPSParam[i].getElementsByTagName('yPosition')[0].childNodes[0].data)
            tmpGPSP.zPosition = float(_GPSParam[i].getElementsByTagName('zPosition')[0].childNodes[0].data)
            tmpGPSP.xVelocity = float(_GPSParam[i].getElementsByTagName('xVelocity')[0].childNodes[0].data)
            tmpGPSP.yVelocity = float(_GPSParam[i].getElementsByTagName('yVelocity')[0].childNodes[0].data)
            tmpGPSP.zVelocity = float(_GPSParam[i].getElementsByTagName('zVelocity')[0].childNodes[0].data)
            tmpGPSP.ZeroTime = _op_.time_dif(self.Imageinfo.imagingTime_start,tmpGPSP.TimeStamp)
            self.GPSParam.append(tmpGPSP)
            
        _ATTIParam = root.getElementsByTagName('ATTIParam')
        for i in range(0,_ATTIParam.__len__()):
            tmpATTP = ATTIParam()
            tmpATTP.TimeStamp = _ATTIParam[i].getElementsByTagName('TimeStamp')[0].childNodes[0].data
            tmpATTP.yawAngle = float(_ATTIParam[i].getElementsByTagName('yawAngle')[0].childNodes[0].data)
            tmpATTP.rollAngle = float(_ATTIParam[i].getElementsByTagName('rollAngle')[0].childNodes[0].data)
            tmpATTP.pitchAngle = float(_ATTIParam[i].getElementsByTagName('pitchAngle')[0].childNodes[0].data)
            self.ATTIParam.append(tmpATTP)
            
        self.Productinfo.NominalResolution = float(root.getElementsByTagName('NominalResolution')[0].childNodes[0].data)
        self.Productinfo.WidthInMeters = float(root.getElementsByTagName('WidthInMeters')[0].childNodes[0].data)
        self.Productinfo.productLevel = int(root.getElementsByTagName('productLevel')[0].childNodes[0].data)
        self.Productinfo.productType = root.getElementsByTagName('productType')[0].childNodes[0].data
        self.Productinfo.productFormat = root.getElementsByTagName('productFormat')[0].childNodes[0].data
        self.Productinfo.productGentime = root.getElementsByTagName('productGentime')[0].childNodes[0].data
        self.Productinfo.productPolar = root.getElementsByTagName('productPolar')[0].childNodes[0].data
        
        self.Processinfo.EphemerisData = root.getElementsByTagName('processinfo')[0].getElementsByTagName('EphemerisData')[0].childNodes[0].data
        self.Processinfo.AttitudeData = root.getElementsByTagName('processinfo')[0].getElementsByTagName('AttitudeData')[0].childNodes[0].data
        self.Processinfo.algorithm = root.getElementsByTagName('processinfo')[0].getElementsByTagName('algorithm')[0].childNodes[0].data
        _Calibration = root.getElementsByTagName('processinfo')[0].getElementsByTagName('CalibrationConst')[0]
        self.Processinfo.Calibration = [_Calibration.getElementsByTagName('HH')[0].childNodes[0].data,\
                                       _Calibration.getElementsByTagName('HV')[0].childNodes[0].data,\
                                       _Calibration.getElementsByTagName('VH')[0].childNodes[0].data,\
                                       _Calibration.getElementsByTagName('VV')[0].childNodes[0].data]
        for i in range(0,self.Processinfo.Calibration.__len__()):
            if self.Processinfo.Calibration[i] == 'NULL':
                self.Processinfo.Calibration[i] = 0
            else:
                self.Processinfo.Calibration[i] = float(self.Processinfo.Calibration[i])
        self.Processinfo.AzFdc0 = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('AzFdc0')[0].childNodes[0].data)
        self.Processinfo.AzFdc1 = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('AzFdc1')[0].childNodes[0].data)
        self.Processinfo.MultilookRange = int(root.getElementsByTagName('processinfo')[0].getElementsByTagName('MultilookRange')[0].childNodes[0].data)
        self.Processinfo.MultilookAzimuth = int(root.getElementsByTagName('processinfo')[0].getElementsByTagName('MultilookAzimuth')[0].childNodes[0].data)
        self.Processinfo.DoIQComp = int(root.getElementsByTagName('processinfo')[0].getElementsByTagName('DoIQComp')[0].childNodes[0].data)
        self.Processinfo.DoChaComp = int(root.getElementsByTagName('processinfo')[0].getElementsByTagName('DoChaComp')[0].childNodes[0].data)
        self.Processinfo.RangeWeightType = int(root.getElementsByTagName('processinfo')[0].getElementsByTagName('RangeWeightType')[0].childNodes[0].data)
        _RangeWeightPara = root.getElementsByTagName('processinfo')[0].getElementsByTagName('RangeWeightPara')[0].childNodes[0].data
        self.Processinfo.RangeWeightPara = [float(_RangeWeightPara.split(',')[0]),float(_RangeWeightPara.split(',')[1])]
        self.Processinfo.AzimuthWeightType = int(root.getElementsByTagName('processinfo')[0].getElementsByTagName('AzimuthWeightType')[0].childNodes[0].data)
        _AzimuthWeightPara = root.getElementsByTagName('processinfo')[0].getElementsByTagName('AzimuthWeightPara')[0].childNodes[0].data
        self.Processinfo.AzimuthWeightPara = [float(_AzimuthWeightPara.split(',')[0]),float(_RangeWeightPara.split(',')[1])]
        self.Processinfo.EarthModel = root.getElementsByTagName('processinfo')[0].getElementsByTagName('EarthModel')[0].childNodes[0].data
        self.Processinfo.ProjectModel = root.getElementsByTagName('processinfo')[0].getElementsByTagName('ProjectModel')[0].childNodes[0].data
        self.Processinfo.incidenceAngleNearRange = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('incidenceAngleNearRange')[0].childNodes[0].data)
        self.Processinfo.incidenceAngleFarRange = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('incidenceAngleFarRange')[0].childNodes[0].data)
        self.Processinfo.RangeLookBandWidth = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('RangeLookBandWidth')[0].childNodes[0].data)
        self.Processinfo.AzimuthLookBandWidth = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('AzimuthLookBandWidth')[0].childNodes[0].data)
        self.Processinfo.TotalProcessedAzimuthBandWidth = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('TotalProcessedAzimuthBandWidth')[0].childNodes[0].data)
        self.Processinfo.DopplerParametersReferenceTime = float(root.getElementsByTagName('processinfo')[0].getElementsByTagName('DopplerParametersReferenceTime')[0].childNodes[0].data)
        _DopplerCentroidCoefficients = root.getElementsByTagName('processinfo')[0].getElementsByTagName('DopplerCentroidCoefficients')[0]
        self.Processinfo.DopplerCentroidCoefficients = [float(_DopplerCentroidCoefficients.getElementsByTagName('d0')[0].childNodes[0].data),\
                                       float(_DopplerCentroidCoefficients.getElementsByTagName('d1')[0].childNodes[0].data),\
                                       float(_DopplerCentroidCoefficients.getElementsByTagName('d2')[0].childNodes[0].data),\
                                       float(_DopplerCentroidCoefficients.getElementsByTagName('d3')[0].childNodes[0].data),\
                                       float(_DopplerCentroidCoefficients.getElementsByTagName('d4')[0].childNodes[0].data)]
        _DopplerRateValuesCoefficients = root.getElementsByTagName('processinfo')[0].getElementsByTagName('DopplerRateValuesCoefficients')[0]
        self.Processinfo.DopplerRateValuesCoefficients = [float(_DopplerRateValuesCoefficients.getElementsByTagName('r0')[0].childNodes[0].data),\
                                       float(_DopplerRateValuesCoefficients.getElementsByTagName('r1')[0].childNodes[0].data),\
                                       float(_DopplerRateValuesCoefficients.getElementsByTagName('r2')[0].childNodes[0].data),\
                                       float(_DopplerRateValuesCoefficients.getElementsByTagName('r3')[0].childNodes[0].data),\
                                       float(_DopplerRateValuesCoefficients.getElementsByTagName('r4')[0].childNodes[0].data)]
        return 0
    def read_data(self,filename):
        dataset = gdal.Open(filename)
        self.Dataset = dataset
        del dataset
        return 0
    def read_incidence(self,incidence_XML):
        incidence_dom = xml.dom.minidom.parse(incidence_XML)
        root = incidence_dom.documentElement

        self.Incidence.nums = int(root.getElementsByTagName('numberofIncidenceValue')[0].childNodes[0].data)
        self.Incidence.stepSize = int(root.getElementsByTagName('stepSize')[0].childNodes[0].data)
        _Value = root.getElementsByTagName('incidenceValue')
        for i in range(0,_Value.__len__()):
            self.Incidence.Value.append(float(_Value[i].childNodes[0].data))
        return 0
    def get_satimaging(self):
        '''
        获取图像成像对应的卫星状态（包括成像前后两个状态与处于成像期间的状态）
        '''
        tmp1 = []
        tmp2 = []
        for i in self.GPSParam:
            tmp1.append(i.ZeroTime)
            if i.ZeroTime >= 0:
                tmp2.append(i.ZeroTime)
        tmp2.append(self.Imageinfo.TimeTotal)
        tmp2.sort()
        tmp3 = tmp2[0:tmp2.index(self.Imageinfo.TimeTotal)+2]
        tmp3.remove(self.Imageinfo.TimeTotal) 
        ind = tmp1.index(tmp3[0])
        tmp3.append(tmp1[ind-1])
        tmp3.sort()
        tmp4=[]
        for j in range(ind-1,ind-1+tmp3.__len__()):
            tmp4.append([tmp3[j-ind+1],self.GPSParam[j].xPosition,self.GPSParam[j].yPosition,self.GPSParam[j].zPosition,self.GPSParam[j].xVelocity,self.GPSParam[j].yVelocity,self.GPSParam[j].zVelocity])
        return tmp4
    def get_satbyT(self,zTime):
        '''
        给定成像时间（以图像第一行为0时间），计算对应卫星状态
        '''
        sati = self.get_satimaging()
        tmp = []
        for i in sati:
            tmp.append(i[0])
        tmp.append(zTime)
        tmp.sort()
        ind = tmp.index(zTime)
        
        t1 = sati[ind-1][0]-sati[0][0]
        t2 = sati[ind][0]-sati[0][0]
        t = zTime-sati[0][0]
        
        sx1 = sati[ind-1][1]
        vx1 = sati[ind-1][4]
        sx2 = sati[ind][1]
        vx2 = sati[ind][4]
        
        sy1 = sati[ind-1][2]
        vy1 = sati[ind-1][5]
        sy2 = sati[ind][2]
        vy2 = sati[ind][5]
        
        sz1 = sati[ind-1][3]
        vz1 = sati[ind-1][6]
        sz2 = sati[ind][3]
        vz2 = sati[ind][6]
        
        A = np.array([[1,t1,t1**2,t1**3],[0,1,2*t1,3*t1**2],[1,t2,t2**2,t2**3],[0,1,2*t2,3*t2**2]])
        Bx = np.transpose(np.array([[sx1,vx1,sx2,vx2]]))
        By = np.transpose(np.array([[sy1,vy1,sy2,vy2]]))
        Bz = np.transpose(np.array([[sz1,vz1,sz2,vz2]]))
        ax = np.linalg.solve(A,Bx).reshape(-1).tolist()
        ay = np.linalg.solve(A,By).reshape(-1).tolist()
        az = np.linalg.solve(A,Bz).reshape(-1).tolist()
        
        Xs = ax[0]+ax[1]*t+ax[2]*t*t+ax[3]*t*t*t
        Xv = ax[1]+2*ax[2]*t+3*ax[3]*t*t
        Xa = 2*ax[2]+6*ax[3]*t
        
        Ys = ay[0]+ay[1]*t+ay[2]*t*t+ay[3]*t*t*t
        Yv = ay[1]+2*ay[2]*t+3*ay[3]*t*t
        Ya = 2*ay[2]+6*ay[3]*t
        
        Zs = az[0]+az[1]*t+az[2]*t*t+az[3]*t*t*t
        Zv = az[1]+2*az[2]*t+3*az[3]*t*t
        Za = 2*az[2]+6*az[3]*t
          
        return [Xs,Ys,Zs,Xv,Yv,Zv,Xa,Ya,Za]
    def get_satAllRow(self):
        '''
        计算图像每一行对应的卫星状态
        '''
        for i in range(0,self.Imageinfo.height):
            zTime = i*self.Imageinfo.TimeInterval
            self.SatEveryRow.append(self.get_satbyT(zTime))
        return 0
    def get_Calib(self,meta_XML):
        '''
        对GF3影像进行辐射定标
        :param meta_XML:元数据文件
        '''
        path = os.path.dirname(meta_XML)
        if self.Productinfo.productLevel != 1:
            print('Do Not Need Calibration')
            return -1
        QualifyValue = {'HH':self.Imageinfo.QualifyValue[0],'HV':self.Imageinfo.QualifyValue[1],'VH':self.Imageinfo.QualifyValue[2],'VV':self.Imageinfo.QualifyValue[3]}
        CalibValue = {'HH':self.Processinfo.Calibration[0],'HV':self.Processinfo.Calibration[1],'VH':self.Processinfo.Calibration[2],'VV':self.Processinfo.Calibration[3]}
        for i in range(self.polarNums):
            time1 = time.time()
            filematch = path + '\\*_' + self.polarMode[i] + '_*[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].tiff'
            filename = glob.glob(filematch)
            Qv = QualifyValue[self.polarMode[i]]
            Cv = CalibValue[self.polarMode[i]]
            if filename.__len__() == 1:
                filein = filename[0]
                if str.find(filein,'L1A') > -1:
                    self.calib_L1A(filein,Qv,Cv)
                time2 = time.time()
                print('Calib ',i,':',time2-time1)
            else:
                print('File Not Found!')
    def calib_L1A(self,filein,QualifyValue,CalibValue):
        '''
        对L1A影像进行辐射定标
        '''
        gdal.SetCacheMax(2**30)
        dataset = gdal.Open(filein)
        data = dataset.ReadAsArray(0,0,dataset.RasterXSize,dataset.RasterYSize,buf_type='float32')
        description = dataset.GetDriver().GetDescription()
        geotrans = dataset.GetGeoTransform()
        projection = dataset.GetProjection()
        dcal = data[0]**2+data[1]**2
        del data
               
        Qtmp = (QualifyValue/32767)**2
        dcal *= Qtmp
        dcal = np.log10(dcal)
        dcal *= 10
        dcal -= CalibValue
        dcal[np.isnan(dcal)] = 0
        dcal[np.isinf(dcal)] = 0
            
        fileout = filein[0:filein.__len__()-5]+'_Calib.tiff'
        _op_.writeTiff(dcal,fileout,description,geotrans,projection)
        del dataset
        del dcal
        return 0    
    
    def help(self):
        print('________________________________________________')
        print('********************GF3_meta********************')
        print('Sensor: lamda,RadarCenterFrequency') 
        print('Wave: centerLookAngle,prf,sampleRate,')
        print(' sampleDelay,groundVelocity') 
        print('lookDirection') 
        print('antennaMode') 
        print('polarMode') 
        print('polarNums') 
        print('Platform: CenterTime,Rs,satVelocity,RollAngle,')
        print(' PitchAngle,YawAngle,Xs,Ys,Zs,Vxs,Vys,Vzs') 
        print('GPSParam: TimeStamp,xPosition,yPosition,')
        print(' zPosition,xVelocity,yVelocity,zVelocity')
        print('ATTIParam: TimeStamp,yawAngle,rollAngle,pitchAngle') 
        print('Productinfo: NominalResolution,WidthInMeters,')
        print(' productLevel,productType,productFormat,')
        print(' productGentime,productPolar') 
        print('Imageinfo:imagingTime_start,imagingTime_end,')
        print(' nearRange,refRange,eqvFs,eqvPRF,center,topLeft,')
        print(' topRight,bottomLeft,bottomRight,width,height,')
        print(' widthspace,heightspace,sceneShift,imagebit,')
        print(' QualifyValue,echoSaturation,TimeTotal,TimeInterval')
        print('Processinfo: EphemerisData,AttitudeData,algorithm,')
        print(' Calibration,AzFdc0,AzFdc1,MultilookRange,')
        print(' MultilookAzimuth,DoIQComp,DoChaComp,RangeWeightType,')
        print(' RangeWeightPara,AzimuthWeightType,AzimuthWeightPara,')
        print(' EarthModel,ProjectModel,incidenceAngleNearRange,')
        print(' incidenceAngleFarRange,RangeLookBandWidth,')
        print(' AzimuthLookBandWidth,TotalProcessedAzimuthBandWidth,')
        print(' DopplerParametersReferenceTime,')
        print(' DopplerCentroidCoefficients,')
        print(' DopplerRateValuesCoefficients') 
        print('SatEveryRow')
        print('************************************************')
        print('________________________________________________')
        print('*****************GF3_incidence******************')
        print('Incidence: nums,stepSize,Value') 
        print('************************************************')
        return 0
        
