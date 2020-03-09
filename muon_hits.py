#!/usr/bin/env python2.7
# -*- coding: utf8 -*-

import time
import datetime
import subprocess
import cStringIO as cstr
import scipy.optimize as optim
import scipy.stats as stats
import numpy as np
import logging as lg
import argparse
import fnmatch
import os

lg.basicConfig(filename='sci.log', level=lg.DEBUG,
               format='%(asctime)s %(message)s')

def parabola(x,k,m,c):
  return k-c*np.square(x-m)

def fits(data,bins,M,N):
  log_data=np.ma.log10(data)
  log_data.filled(0)
  amax=np.argmax(log_data,axis=2)
  nfits=0
  pdata=np.zeros([2,M,N])
  for k in range(0,M):
    for j in range(0,N):
      xdata=np.arange(amax[k,j]-4,amax[k,j]+5+1)
      ydata=log_data[k,j,xdata]
      if np.all(ydata):
        try:
          popt,pcov=optim.curve_fit(parabola,xdata,ydata)
          perr=np.sqrt(np.diag(pcov))
          if np.all(perr<1.0):
            m_ADC=np.uint16(np.rint(popt[1]))
            pdata[0,k,j]=m_ADC
            pdata[1,k,j]=np.sum(data[k,j,m_ADC+290:])
            nfits+=1
          else:
            m_ADC=amax[k,j]
            pdata[0,k,j]=m_ADC
            pdata[1,k,j]=np.sum(data[k,j,m_ADC+290:])
        except RuntimeError:
          m_ADC=amax[k,j]
          pdata[0,k,j]=m_ADC
          pdata[1,k,j]=np.sum(data[k,j,m_ADC+290:])
  return pdata,nfits

def count_rate(home,nfile,adc,ch_num):
  fs=[]
  for yx in range(27,29):
    name='{0}/{1}.{2}.gz'.format(home,nfile[:-6],yx)
    p=subprocess.Popen(['/usr/bin/zcat',name],stdout=subprocess.PIPE)
    f=io_method(p.communicate()[0])
    assert p.returncode==0
    f.readline()
    fs.append(f)
  line0=fs[0].readline()
  line1=fs[1].readline()
  M=len(line0.split())+len(line1.split())
  eventnum=0
  while M==548:
    line_test0=fs[0].readline()
    line_test1=fs[1].readline()
    k0=np.fromstring(line0,dtype=np.uint16,count=265,sep=' ')[9:].reshape(64,4)
    k1=np.fromstring(line1,dtype=np.uint16,count=265,sep=' ')[9:].reshape(64,4)
    adc[0,0,ch_num,k0[:,0]]+=1.0
    adc[0,1,ch_num,k0[:,1]]+=1.0
    adc[0,2,ch_num,k0[:,2]]+=1.0
    adc[0,3,ch_num,k0[:,3]]+=1.0
    adc[1,0,ch_num,k1[:,0]]+=1.0
    adc[1,1,ch_num,k1[:,1]]+=1.0
    adc[1,2,ch_num,k1[:,2]]+=1.0
    adc[1,3,ch_num,k1[:,3]]+=1.0
    M=len(line_test0.split())+len(line_test1.split())
    line0=line_test0
    line1=line_test1
    eventnum+=1
  f.close()
  adc_top=adc[0,:,:,:]
  adc_bot=adc[1,:,:,:]
  return adc_top,adc_bot,eventnum

parser=argparse.ArgumentParser()
parser.add_argument('home', help='home',type=str)
parser.add_argument('year', help='year number',type=int)
parser.add_argument('month', help='month number',type=int)
parser.add_argument('bdate', help='begin date',type=int)
parser.add_argument('edate', help='end date',type=int)

args=parser.parse_args()
home=args.home
year=args.year
month=args.month
bdate=args.bdate
edate=args.edate
d0=datetime.datetime(2000+year,month,bdate)

# for python3 use import io instead of cStringIO
# and io.BytesIO instead cstr.StringIO
io_method=cstr.StringIO
lg.info('begin date: {0}, end date: {1}'.format(bdate,edate))

t0=time.time()
time_stamp=[]
wdates=np.arange(bdate,edate+1)
ndays=np.size(wdates)
nhoras=24
muon_files=os.listdir(home)
for d in wdates:
  mpat='SN.0.{0}{1:02d}{2:02d}*.27.gz'.format(year,month,d)
  for mfile in muon_files:
    if fnmatch.fnmatch(mfile,mpat):
      time_stamp.append(mfile)

Nch=64
Mbins=2**12
FEBnum=4
ebins=np.arange(0,Mbins)
Nchans=np.arange(0,Nch)
hit_data=np.zeros([ndays,nhoras,9])
for ts in time_stamp:
  cday=int(ts[9:11])
  day_num=cday-bdate
  chour=np.int(ts[11:13])
  adc=np.zeros([2,FEBnum,Nch,Mbins])
  lg.info('Procesando archivo: {0}'.format(ts))
  adc_top,adc_bot,eventnum=count_rate(home,ts,adc,Nchans)
  ped_m0,nfits0=fits(adc_top,ebins,FEBnum,Nch)
  ped_m1,nfits1=fits(adc_bot,ebins,FEBnum,Nch)
  hit_data[day_num,chour,0]+=eventnum
  hit_data[day_num,chour,1:5]+=np.sum(ped_m0[1,:,:],axis=1)
  hit_data[day_num,chour,5:9]+=np.sum(ped_m1[1,:,:],axis=1)

sep=' '
f=open('hits-muons-nov.txt','w')
for m in range(0,ndays):
  for j in range(0,nhoras):
    dt=datetime.timedelta(days=m,hours=j)
    date=(d0+dt).strftime('%Y-%m-%d %H:00 ')
    dataA='{0} {1} '.format(nfits0,nfits1)
    dataH=sep.join(str(x) for x in np.uint32(hit_data[m,j,:]))
    f.write(date+dataA+dataH)
    f.write('\n')
f.close()

t1=time.time()
dt=str(datetime.timedelta(seconds=(t1-t0)))[0:7]
lg.info('Tiempo total {0} hrs.'.format(dt))
