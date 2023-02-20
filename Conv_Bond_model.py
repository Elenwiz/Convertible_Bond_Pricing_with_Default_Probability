# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:10:14 2022

@author: Chenyue
"""
import numpy as np
import time
from datetime import datetime, timedelta
from chinese_calendar import is_workday, is_holiday
import tushare as ts
from scipy.interpolate import interp1d

YearDays = 250
# input your tushare token
pro = ts.pro_api()


def ConvertibleBond(S0, K, T, r, spread, sigma, UnLockDate, PayDates, CouponRates, Nodes=1000, dividend=0, Lambda=0, Theta=0, R=0):
    '''
    S0: Current Stock Price
    K: Conversion Price
    T: Current Time
    r: risk-free interst rate
    spread: credit spread
    dividend: the dividend ratio of the underlying stock(0 by default)
    UnLockDate: Date when the conversion starts to be allowed
    PayDates: a list containing dates that the bond pays out coupon cash flow
    CouponRates: a list containing coupon amounts in every paydate
    Lambda: the Poisson process's parameter, which describes the default probability (0 by default)
    Theta: the percentage the stock price drops immediately after the default (0 by default)
    '''

    ResidualDays = (PayDates[-1] - T).days
    # If residual days are not that much, we choose to shorten the time step
    NumInDay = int(Nodes/ResidualDays) + 1
    TreeNum = NumInDay * ResidualDays
    DaysToPay = list(map(lambda x: (x - T).days * NumInDay, PayDates))
    DaysToUnLock = (UnLockDate - T).days * NumInDay

    CRatio = 100/K
    dt = 1/365/NumInDay
    b = r - dividend


    u = np.exp(sigma * np.sqrt(dt))
    d = 1/u
    p_0 = 1 - np.exp(-Lambda * dt)
    p_u = (np.exp(b * dt) - np.exp(-Lambda * dt) * d - (1 - Theta) * p_0)/(u - d)


    mu = np.arange(TreeNum + 1)
    mu = np.resize(mu, (TreeNum + 1, TreeNum + 1))
    md = np.transpose(mu)
    mu = u ** (mu - md)
    md = d ** md
    S = S0 * mu * md


    #Valuation by Backwards Induction
    BondFinalValue = 100*(1+CouponRates[-1])
    ConvertValue = S * CRatio
    DConvertValue = ConvertValue * (1-Theta)
    DefaultValue = np.where(DConvertValue > R * 100, DConvertValue, R*100)
    DConvertProb = np.where(DefaultValue == DConvertValue, 1, 0)
    

    CBondValue = np.zeros((TreeNum+1, TreeNum+1))  # value of CB
    HoldValue = np.zeros((TreeNum+1, TreeNum+1))  # H_t
    ConvertProb = np.zeros((TreeNum+1, TreeNum+1))  # p_t
    DfMatrix = np.zeros((TreeNum+1, TreeNum+1))  # discount factor
    
    
    CBondValue[:,-1] = (1-p_0) * np.where(ConvertValue[:,-1] >= BondFinalValue, ConvertValue[:,-1], BondFinalValue) + p_0 * DefaultValue[:,-1]
    HoldValue[:,-1] = CBondValue[:,-1]
    ConvertProb[:,-1] = np.where(ConvertValue[:,-1] >= BondFinalValue, 1, 0)
    DfMatrix[:,-1] = ConvertProb[:,-1]*np.exp(-b*dt) + (1-ConvertProb[:,-1])*np.exp(-(b+spread)*dt)
    
    
    cnt = 0
    for i in range(TreeNum-1,-1,-1):
        HoldValue[0:TreeNum-cnt,i] = p_u*CBondValue[0:TreeNum-cnt,i+1]*DfMatrix[0:TreeNum-cnt,i+1] + (1-p_0-p_u)*CBondValue[1:TreeNum-cnt+1,i+1]*DfMatrix[1:TreeNum-cnt+1,i+1] + p_0*DefaultValue[0:TreeNum-cnt,i]*DfMatrix[0:TreeNum-cnt,i+1]
        # HoldValue = np.round(HoldValue, 10)
        ConvertProb[0:TreeNum-cnt,i] = p_u*ConvertProb[0:TreeNum-cnt,i+1] + (1-p_0-p_u)*ConvertProb[1:TreeNum-cnt+1,i+1] + p_0 * DConvertProb[0:TreeNum-cnt,i]
        TempCoupon = 100*CouponRates[np.nonzero(np.array(DaysToPay) == i)[0][0]] if i in DaysToPay else 0
        if i < DaysToUnLock:
            CBondValue[0:TreeNum-cnt,i] = (1 - p_0) * (HoldValue[0:TreeNum-cnt,i] + TempCoupon) + p_0 * 100 * R
        else:
            CBondValue[0:TreeNum-cnt,i] = (1 - p_0) * (np.where(HoldValue[0:TreeNum-cnt,i] + TempCoupon > ConvertValue[0:TreeNum-cnt,i],HoldValue[0:TreeNum-cnt,i]+TempCoupon,ConvertValue[0:TreeNum-cnt,i])) + p_0 * DefaultValue[0:TreeNum-cnt,i]
            ConvertProb[0:TreeNum-cnt,i] = (1 - p_0) * (np.where(HoldValue[0:TreeNum-cnt,i] + TempCoupon > ConvertValue[0:TreeNum-cnt,i],ConvertProb[0:TreeNum-cnt,i],1)) + p_0
        DfMatrix[0:TreeNum-cnt,i] = ConvertProb[0:TreeNum-cnt,i]*np.exp(-b*dt) + (1-ConvertProb[0:TreeNum-cnt,i])*np.exp(-(b+spread)*dt)
        cnt += 1

    return CBondValue[0, 0]


def CBProb(S0, K, T, r, spread, sigma, UnLockDate, PayDates, CouponRates, dividend=0, Lambda=0, Theta=0):
    '''
    S0: Current Stock Price
    K: Conversion Price
    T: Current Time
    r: risk-free interst rate
    spread: credit spread
    dividend: the dividend ratio of the underlying stock(0 by default)
    UnLockDate: Date when the conversion starts to be allowed
    PayDates: a list containing dates that the bond pays out coupon cash flow
    CouponRates: a list containing coupon amounts in every paydate
    Lambda: the Poisson process's parameter, which describes the default probability (0 by default)
    Theta: the percentage the stock price drops immediately after the default (0 by default)
    '''

    ResidualDays = (PayDates[-1] - T).days
    # If residual days are not that much, we choose to shorten the time step
    NumInDay = int(1000/ResidualDays) + 1
    TreeNum = ResidualDays * NumInDay
    DaysToPay = list(map(lambda x: (x - T).days * NumInDay, PayDates))
    DaysToUnLock = (UnLockDate - T).days * NumInDay

    CRatio = 100/K
    dt = 1/365/NumInDay
    b = r - dividend


    u = np.exp(sigma * np.sqrt(dt))
    d = 1/u
    p_0 = 1 - np.exp(-Lambda * dt)
    p_u = (np.exp(b * dt) - np.exp(-Lambda * dt) * d - (1 - Theta) * p_0)/(u - d)


    mu = np.arange(TreeNum + 1)
    mu = np.resize(mu, (TreeNum + 1, TreeNum + 1))
    md = np.transpose(mu)
    mu = u ** (mu - md)
    md = d ** md
    S = S0 * mu * md


    #Valuation by Backwards Induction
    BondFinalValue = 100*(1+CouponRates[-1])
    ConvertValue = S * CRatio
    DConvertValue = ConvertValue * (1-Theta)


    CBondValue = np.zeros((TreeNum+1, TreeNum+1))  # value of CB
    HoldValue = np.zeros((TreeNum+1, TreeNum+1))  # H_t
    ConvertProb = np.zeros((TreeNum+1, TreeNum+1))  # p_t
    DfMatrix = np.zeros((TreeNum+1, TreeNum+1))  # discount factor
    
    
    CBondValue[:,-1] = np.where(ConvertValue[:,-1] >= BondFinalValue, ConvertValue[:,-1], BondFinalValue)
    HoldValue[:,-1] = CBondValue[:,-1]
    ConvertProb[:,-1] = np.where(ConvertValue[:,-1] >= BondFinalValue, 1, 0)
    DfMatrix[:,-1] = ConvertProb[:,-1]*np.exp(-b*dt) + (1-ConvertProb[:,-1])*np.exp(-(b+spread)*dt)
    
    
    cnt = 0
    for i in range(TreeNum-1,-1,-1):
        HoldValue[0:TreeNum-cnt,i] = p_u*CBondValue[0:TreeNum-cnt,i+1]*DfMatrix[0:TreeNum-cnt,i+1] + (np.exp(-Lambda * dt)-p_u)*CBondValue[1:TreeNum-cnt+1,i+1]*DfMatrix[1:TreeNum-cnt+1,i+1] + p_0*DConvertValue[0:TreeNum-cnt,i]*DfMatrix[0:TreeNum-cnt,i+1]
        # HoldValue = np.round(HoldValue, 10)
        ConvertProb[0:TreeNum-cnt,i] = p_u*ConvertProb[0:TreeNum-cnt,i+1] + (1-p_0-p_u)*ConvertProb[1:TreeNum-cnt+1,i+1] + p_0
        TempCoupon = 100*CouponRates[np.nonzero(np.array(DaysToPay) == i)[0][0]] if i in DaysToPay else 0
        if i < DaysToUnLock:
            CBondValue[0:TreeNum-cnt,i] = (1 - p_0) * (HoldValue[0:TreeNum-cnt,i] + TempCoupon) + p_0 * DConvertValue[0:TreeNum-cnt,i]
        else:
            CBondValue[0:TreeNum-cnt,i] = (1 - p_0) * (np.where(HoldValue[0:TreeNum-cnt,i] + TempCoupon > ConvertValue[0:TreeNum-cnt,i],HoldValue[0:TreeNum-cnt,i]+TempCoupon,ConvertValue[0:TreeNum-cnt,i])) + p_0 * DConvertValue[0:TreeNum-cnt,i]
            ConvertProb[0:TreeNum-cnt,i] = (1 - p_0) * (np.where(HoldValue[0:TreeNum-cnt,i] + TempCoupon > ConvertValue[0:TreeNum-cnt,i],ConvertProb[0:TreeNum-cnt,i],1)) + p_0
        DfMatrix[0:TreeNum-cnt,i] = ConvertProb[0:TreeNum-cnt,i]*np.exp(-b*dt) + (1-ConvertProb[0:TreeNum-cnt,i])*np.exp(-(r+spread)*dt)
        cnt += 1

    return ConvertProb[0, 0]


def tradedays(start,end):
    '''
    start: start time
    end: end time
    '''
    # deal with the date format
    if type(start) == str:
        start = datetime.strptime(start,'%Y-%m-%d').date()
    if type(end) == str:
        end = datetime.strptime(end,'%Y-%m-%d').date()
    if start > end:
        start,end = end,start
        
    counts = 0
    while True:
        if start > end:
            break
        try:
            if is_holiday(start) or start.weekday()==5 or start.weekday()==6:
                start += timedelta(days=1)
                continue
            counts += 1
            start += timedelta(days=1)
        except:
            counts = end - start
            counts = counts.days
            break
    return counts


def CalVol(uncode,VDate,EndDate):
    TDateCount = np.minimum(np.maximum(tradedays(VDate, EndDate),20),250)

    if type(VDate) == str:
        VDate = datetime.strptime(VDate,'%Y-%m-%d')
    if type(EndDate) == str:
        EndDate = datetime.strptime(EndDate,'%Y-%m-%d')


    delta1 = timedelta(days = int(TDateCount))
    CalBeginDate = VDate - delta1
    
    while tradedays(CalBeginDate, VDate) < TDateCount:
        CalBeginDate = CalBeginDate - timedelta(days = 1)

    CalBeginDate = CalBeginDate.strftime('%Y%m%d')

    TradeSeri = pro.daily(ts_code=uncode, start_date=CalBeginDate, end_date=VDate.strftime('%Y%m%d'), fields='close')
    PSeri = TradeSeri['close']
    LogReturn = np.log(PSeri[1:])-np.log(PSeri.shift()[1:])
    return round(np.std(LogReturn)*np.sqrt(YearDays),4)


def Calyieldrate(Table, rating, date, duration):
    '''
    Table: DataFrame used to calculate the yield rate, the index of which is date
    rating: string, e.g 'AAA'
    date: calculate date, e.g '2022-01-02'
    duration: the rest duration (in year)
    '''
    y_col = []
    title = rating + ':'
    # print(title)
    for i in Table.columns:
        if title in i[:len(title)]:
            y_col.append(i)
    # print(y_col)
    x = [0, 0.083, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]
    y = Table.loc[date, :][y_col]
    # print(y)
    interp_func = interp1d(x, y)
    return interp_func(duration)


def Sum(a):
    count = []
    for i in range(1, len(a)+1):
        count.append(sum(a[:i]))
    return count


def GetConvertPrice(ts_code, trade_date):
    price_chg = pro.cb_price_chg(ts_code=ts_code)
    x = Sum([x>=trade_date for x in price_chg.change_date])
    if sum(x) == 0 or sum(x) == 1:
        convertprice = price_chg['convertprice_aft'][len(x)-1]
    else:
        Index = x.index(1) - 1
        convertprice = price_chg['convertprice_aft'][Index]
    if convertprice is None or np.isnan(convertprice):
        convertprice = price_chg['convert_price_initial'][0]
    return convertprice


if __name__=='__main__':


    uncode = '001979.SZ'
    VDate = datetime(2020,7,28)
    EndDate = datetime(2022,3,13)
    sigma = CalVol(uncode, VDate, EndDate)
    print(sigma)


    S0 = 17.08
    K = 22.08
    T = datetime(2020,7,28)
    rf = 0.02807
    spread = 0.002227
    UnLockDate = datetime(2019,9,16)
    PayDates = [datetime(2020,3,13),datetime(2021,3,13),datetime(2022,3,13)]
    CouponRate = [0.001,0.001,0.061]

    P = []
    for i in S0:
        start = time.time()
        price = ConvertibleBond(i,K,T,rf,spread,sigma,UnLockDate,PayDates,CouponRate,Lambda=(1-np.exp(-spread)),Theta=0.7)
        end = time.time()
        print(end - start)
        print(price)
        P.append(price)
    print(P)
